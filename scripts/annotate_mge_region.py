#!/usr/bin/env python3
"""BLAST each paired tRNA against the genomic window next to it.

The raw hits are the attL candidates; ``filter_confident_ie.py`` picks attL
from them.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import tempfile

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from genome import contig_record, read_genome
from logger import get_logger
from v3ps_filters import DEFAULT_ATTL_WINDOW_BP, load_thresholds

logger = get_logger("annotate_mge_region")

BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]
"""Column order for BLAST tabular output (format 6)."""


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
    query_fa: str,
    out_tsv: str,
    tmp_dir: str = ".",
    window_size: int = DEFAULT_ATTL_WINDOW_BP,
) -> None:
    """BLAST every tRNA query in ``query_fa`` against its window and write all hits.

    The window is ``window_size`` bp on the 3' side of the tRNA. Each hit row
    carries the query id, ``integrase_id``, ``contig`` and ``wstart``, the
    1-based contig coordinate of the window's first base.

    Args:
        genome_fasta: Path to the sample genome FASTA file.
        query_fa: Path to ``mge_query.fa`` (tRNA sequences used as BLAST queries).
        out_tsv: Output path for the raw hits (``mge_blast_raw.tsv``).
        tmp_dir: Directory for temporary BLAST files.
        window_size: Length of the searched window in bp.
    """
    genome = read_genome(genome_fasta)
    frames: list[pd.DataFrame] = []

    for rec in SeqIO.parse(query_fa, "fasta"):
        try:
            integrase_id, contig, coord_range, strand = rec.id.strip().split(":")
            trna_start, trna_end = map(int, coord_range.split("-"))
        except ValueError:
            logger.error(f"Malformed FASTA header: {rec.id}")
            continue

        if strand == "+":
            window_start, window_end = trna_end, trna_end + window_size
        else:
            window_start, window_end = max(1, trna_start - window_size), trna_start
        contig_seq = contig_record(genome, contig).seq
        subject = contig_seq[window_start - 1 : min(len(contig_seq), window_end)]

        hits = run_blast_on_region(
            SeqRecord(rec.seq, id=rec.id, description=""),
            subject,
            tmp_dir,
        )
        if hits.empty:
            logger.info(f"No BLAST hits for {rec.id}")
            continue
        hits["integrase_id"] = integrase_id
        hits["contig"] = contig
        hits["wstart"] = window_start
        frames.append(hits)

    if frames:
        pd.concat(frames).to_csv(out_tsv, sep="\t", index=False)
        logger.info(f"{sum(len(df) for df in frames)} raw hits saved to {out_tsv}")
    else:
        open(out_tsv, "w").close()
        logger.info("No raw BLAST hits collected.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="BLAST tRNA queries against the genomic window next to each tRNA."
    )
    parser.add_argument("--ffn", required=True, help="Genome FASTA path.")
    parser.add_argument("--query", required=True, help="mge_query.fa path.")
    parser.add_argument("--out_tsv", required=True, help="Output mge_blast_raw.tsv path.")
    parser.add_argument("--tmp_dir", default=".", help="Temporary directory for BLAST.")
    parser.add_argument(
        "--params", "--config", dest="params", default="search_params.yaml",
        help="Search parameters (search_params.yaml).",
    )
    cli_args = parser.parse_args()
    thresholds = load_thresholds(cli_args.params)
    main(
        cli_args.ffn, cli_args.query, cli_args.out_tsv, cli_args.tmp_dir,
        window_size=thresholds.attl_window_bp,
    )
