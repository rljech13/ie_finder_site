#!/usr/bin/env python3
"""Run local BLAST to identify attachment-site candidates near integrase-tRNA pairs."""

from __future__ import annotations

import argparse
import os
import tempfile

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from logger import Logger
from v3ps_filters import (
    DEFAULT_ATTL_CANDIDATE_MIN_BP,
    DEFAULT_ATTL_SELECT_BY,
    DEFAULT_ATTL_WINDOW_BP,
    DEFAULT_SHIFT,
    load_thresholds,
    select_v3ps_strict_hit,
)

logger = Logger(name="annotate_mge_region", level=Logger.Level.INFO).get_logger()

BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]
"""Column order for BLAST tabular output (format 6)."""


def load_integrases(integrases_file: str) -> pd.DataFrame:
    """Load integrase hit summary table.

    Args:
        integrases_file: Path to ``integrase_hits_summary.tsv``.

    Returns:
        Parsed integrase DataFrame.
    """
    return pd.read_csv(integrases_file, sep="\t")


def get_integrase_record(df: pd.DataFrame, integrase_id: str) -> pd.Series | None:
    """Return the integrase row matching ``integrase_id``.

    Args:
        df: Integrase summary DataFrame.
        integrase_id: Integrase ORF identifier.

    Returns:
        Matching row as a Series, or None when the identifier is absent.
    """
    rec = df[df["orf_id"] == integrase_id]
    if rec.empty:
        logger.error(f"Integrase ID {integrase_id} not found in integrase table")
        return None
    return rec.iloc[0]


def extract_subject_region(
    genome_fasta: str,
    contig: str,
    region_start: int,
    region_end: int,
) -> Seq | None:
    """Extract a genomic subsequence used as the BLAST subject.

    Args:
        genome_fasta: Path to the whole-genome FASTA file.
        contig: Contig identifier.
        region_start: Window start coordinate (1-based inclusive).
        region_end: Window end coordinate (1-based inclusive).

    Returns:
        Subsequence as a Biopython ``Seq``, or None when the contig is not found.
    """
    recs = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    rec = recs.get(contig)
    matched_key = contig
    if not rec:
        for key, candidate in recs.items():
            if contig in key:
                rec = candidate
                matched_key = key
                logger.info(f"Using contig {matched_key} for match {contig}")
                break
    if not rec:
        logger.error(f"Contig {contig} not found")
        return None
    seq_len = len(rec.seq)
    start = max(1, region_start)
    end = min(seq_len, region_end)
    return rec.seq[start - 1 : end]


def run_blast_on_region(
    query_rec: SeqRecord,
    subject_seq: Seq,
    tmp_dir: str = ".",
) -> pd.DataFrame:
    """Run nucleotide BLAST of one tRNA query against a genomic window.

    Args:
        query_rec: tRNA query sequence record.
        subject_seq: Genomic window sequence.
        tmp_dir: Directory for temporary BLAST database files.

    Returns:
        BLAST hits as a DataFrame with columns listed in ``BLAST_COLS``.
        Returns an empty DataFrame when BLAST produces no alignments.
    """
    import subprocess

    with tempfile.TemporaryDirectory(dir=tmp_dir) as tmp:
        qfa = os.path.join(tmp, "query.fa")
        sfa = os.path.join(tmp, "subject.fa")
        bout = os.path.join(tmp, "blast.tsv")

        SeqIO.write([query_rec], qfa, "fasta")
        SeqIO.write([SeqRecord(subject_seq, id="subject", description="")], sfa, "fasta")

        subprocess.run(
            ["makeblastdb", "-in", sfa, "-dbtype", "nucl"],
            check=True,
            capture_output=True,
            text=True,
        )
        subprocess.run(
            [
                "blastn", "-query", qfa, "-db", sfa,
                "-outfmt", f"6 {' '.join(BLAST_COLS)}",
                "-word_size", "4", "-dust", "no",
                "-out", bout,
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        try:
            return pd.read_csv(bout, sep="\t", header=None, names=BLAST_COLS)
        except pd.errors.EmptyDataError:
            return pd.DataFrame(columns=BLAST_COLS)


def main(
    genome_fasta: str,
    integrases_file: str,
    query_fa: str,
    out_tsv: str,
    tmp_dir: str = ".",
    window_size: int = DEFAULT_ATTL_WINDOW_BP,
    shift: int = DEFAULT_SHIFT,
    min_len_bp: int = DEFAULT_ATTL_CANDIDATE_MIN_BP,
    select_by: str = DEFAULT_ATTL_SELECT_BY,
) -> None:
    """Run attachment-site BLAST for each integrase-tRNA query and write results.

    For each tRNA query in ``query_fa``, BLAST is run against a genomic window
    around the paired integrase. The best 3'-anchored, strand-consistent hit
    of at least ``min_len_bp`` is written to ``out_tsv``: the highest bitscore,
    or the longest hit when ``select_by`` is ``"length"``. All raw hits are
    written to a sibling file with the ``_raw.tsv`` suffix.

    Args:
        genome_fasta: Path to the sample genome FASTA file.
        integrases_file: Path to ``integrase_hits_summary.tsv``.
        query_fa: Path to ``mge_query.fa`` (tRNA sequences used as BLAST queries).
        out_tsv: Output path for filtered best hits (``mge_blast.tsv``).
        tmp_dir: Directory for temporary BLAST files.
        select_by: ``"bitscore"`` or ``"length"``; see ``select_v3ps_strict_hit``.
    """
    integrases = load_integrases(integrases_file)
    filtered_records: list[dict] = []
    raw_frames: list[pd.DataFrame] = []

    for rec in SeqIO.parse(query_fa, "fasta"):
        try:
            integrase_id, contig, coord_range, strand = rec.id.strip().split(":")
            trna_start, trna_end = map(int, coord_range.split("-"))
        except ValueError:
            logger.error(f"Malformed FASTA header: {rec.id}")
            continue

        integrase = get_integrase_record(integrases, integrase_id)
        if integrase is None:
            continue

        if strand == "+":
            window_start, window_end = trna_end, trna_end + window_size
        else:
            window_start, window_end = max(1, trna_start - window_size), trna_start

        subject = extract_subject_region(genome_fasta, contig, window_start, window_end)
        if subject is None:
            logger.error(f"No BLAST window extracted for {rec.id}")
            continue

        blast_df = run_blast_on_region(
            SeqRecord(rec.seq, id=rec.id, description=""),
            subject,
            tmp_dir,
        )
        if blast_df.empty:
            logger.info(f"No BLAST hits for {rec.id}")
            continue

        blast_df["integrase_id"] = integrase_id
        blast_df["contig"] = contig
        blast_df["wstart"] = window_start
        raw_frames.append(blast_df.copy())

        best = select_v3ps_strict_hit(
            blast_df, len(rec.seq), strand, shift, min_len_bp, select_by
        )
        if best is None:
            logger.info(f"No strict attachment-site hits for {rec.id}")
            continue

        filtered_records.append({
            "integrase_id": integrase_id,
            "contig": contig,
            "hit_start": best["attL_abs_lo"],
            "hit_end": best["attL_abs_hi"],
            "pident": best["attL_pident"],
            "length": best["best_attl_len_bp"],
            "mismatch": best["attL_mismatch"],
            "gapopen": best["attL_gapopen"],
            "evalue": best["attL_evalue"],
            "bitscore": best["attL_bitscore"],
            "qstart": best["qstart"],
            "qend": best["qend"],
        })

    if filtered_records:
        pd.DataFrame(filtered_records).to_csv(out_tsv, sep="\t", index=False)
        logger.info(f"{len(filtered_records)} filtered hits saved to {out_tsv}")
    else:
        open(out_tsv, "w").close()
        logger.info("No valid BLAST hits; wrote empty filtered table.")

    raw_path = out_tsv.replace(".tsv", "_raw.tsv")
    if raw_frames:
        pd.concat(raw_frames).to_csv(raw_path, sep="\t", index=False)
        logger.info(f"{sum(len(df) for df in raw_frames)} raw hits saved to {raw_path}")
    else:
        open(raw_path, "w").close()
        logger.info("No raw BLAST hits collected.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="BLAST tRNA queries against genomic windows near integrases."
    )
    parser.add_argument("--ffn", required=True, help="Genome FASTA path.")
    parser.add_argument("--integrases", required=True, help="integrase_hits_summary.tsv path.")
    parser.add_argument("--query", required=True, help="mge_query.fa path.")
    parser.add_argument("--out_tsv", required=True, help="Output mge_blast.tsv path.")
    parser.add_argument("--tmp_dir", default=".", help="Temporary directory for BLAST.")
    parser.add_argument(
        "--params", "--config", dest="params", default="search_params.yaml",
        help="Search parameters (search_params.yaml).",
    )
    cli_args = parser.parse_args()
    thresholds = load_thresholds(cli_args.params)
    main(
        cli_args.ffn, cli_args.integrases, cli_args.query, cli_args.out_tsv, cli_args.tmp_dir,
        window_size=thresholds.attl_window_bp,
        shift=thresholds.shift,
        min_len_bp=thresholds.attl_candidate_min_bp,
        select_by=thresholds.attl_select_by,
    )
