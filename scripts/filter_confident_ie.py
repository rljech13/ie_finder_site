#!/usr/bin/env python3
"""Apply confidence filters to integrative element candidates for one sample."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from logger import Logger
from sequence_qc import parse_fasta_lengths_and_n_flags
from v3ps_filters import (
    FilterThresholds,
    build_qseqid,
    evaluate_ie_candidate,
    load_thresholds,
    parse_orfs_gff,
)

logger = Logger(name="filter_confident_ie").get_logger()


def closest_trna_rows(trna_path: Path) -> pd.DataFrame:
    """Load integrase-tRNA pairs, retaining the closest tRNA per integrase.

    Args:
        trna_path: Path to ``integrase_trna.tsv``.

    Returns:
        DataFrame with at most one row per ``integrase_id``.
    """
    df = pd.read_csv(trna_path, sep="\t", dtype={"contig": str})
    if df.empty:
        return df
    if "distance" in df.columns:
        df = df.loc[df.groupby("integrase_id")["distance"].idxmin()]
    return df


def integrase_coords(hits_path: Path) -> dict[str, tuple[int, int]]:
    """Load integrase nucleotide coordinates keyed by ORF identifier.

    Args:
        hits_path: Path to ``integrase_hits_summary.tsv``.

    Returns:
        Dictionary mapping ``orf_id`` to ``(start, end)`` (1-based inclusive).
    """
    df = pd.read_csv(hits_path, sep="\t")
    out: dict[str, tuple[int, int]] = {}
    for _, row in df.iterrows():
        out[str(row["orf_id"])] = (int(row["start"]), int(row["end"]))
    return out


def load_raw_blast(raw_path: Path) -> pd.DataFrame:
    """Load raw BLAST output, returning an empty frame when the file is missing.

    Args:
        raw_path: Path to ``mge_blast_raw.tsv``.

    Returns:
        Parsed BLAST DataFrame, or an empty DataFrame when the file is absent.
    """
    if not raw_path.is_file() or raw_path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(raw_path, sep="\t", dtype={"contig": str})


def filter_sample(
    *,
    sample: str,
    trna_path: Path,
    integrase_hits_path: Path,
    raw_blast_path: Path,
    mge_fa_path: Path,
    mge_gbk_path: Path,
    orfs_gff_path: Path,
    fasta_path: Path,
    out_fa: Path,
    out_gbk: Path,
    out_audit: Path,
    thresholds: FilterThresholds,
) -> int:
    """Filter integrative element candidates and write confident outputs for one sample.

    Args:
        sample: Sample identifier (genome name).
        trna_path: Path to ``integrase_trna.tsv``.
        integrase_hits_path: Path to ``integrase_hits_summary.tsv``.
        raw_blast_path: Path to ``mge_blast_raw.tsv``.
        mge_fa_path: Path to ``mge_region.fa`` containing candidate elements.
        mge_gbk_path: Path to ``mge_annotated.gbk``.
        orfs_gff_path: Path to Prodigal ``orfs.gff``.
        fasta_path: Path to the assembly FASTA. Exact attL runs are read from it.
        out_fa: Output path for ``ie_confident.fa``.
        out_gbk: Output path for ``ie_confident.gbk``.
        out_audit: Output path for ``ie_filter_audit.tsv``.
        thresholds: Confidence filter thresholds.

    Returns:
        Number of integrase identifiers that passed all confidence filters.
    """
    trna_df = closest_trna_rows(trna_path)
    int_coords = integrase_coords(integrase_hits_path)
    raw_blast = load_raw_blast(raw_blast_path)
    raw_by_q = raw_blast.groupby("qseqid") if not raw_blast.empty else None

    ie_lengths: dict[str, int] = {}
    ie_has_n: dict[str, bool] = {}
    ie_seqs: dict[str, str] = {}
    if mge_fa_path.is_file() and mge_fa_path.stat().st_size > 0:
        ie_lengths, ie_has_n = parse_fasta_lengths_and_n_flags(mge_fa_path)
        for rec in SeqIO.parse(mge_fa_path, "fasta"):
            ie_seqs[rec.id] = str(rec.seq).upper()

    # attL / attR coordinates as written by extract_mge_regions.py
    sites_path = mge_fa_path.parent / "attachment_sites.tsv"
    sites: dict[str, pd.Series] = {}
    if sites_path.is_file() and sites_path.stat().st_size > 0:
        sites_df = pd.read_csv(sites_path, sep="\t")
        if not sites_df.empty:
            sites_df["integrase_id"] = sites_df["integrase_id"].astype(str)
            sites = sites_df.drop_duplicates("integrase_id").set_index("integrase_id").to_dict("index")

    cds_by_contig = parse_orfs_gff(orfs_gff_path)
    contig_seqs: dict[str, str] = {}
    if fasta_path.is_file() and fasta_path.stat().st_size > 0:
        contig_seqs = {
            rec.id: str(rec.seq).upper()
            for rec in SeqIO.parse(fasta_path, "fasta")
        }
    audit_rows: list[dict] = []

    if trna_df.empty:
        out_fa.write_text("")
        out_gbk.write_text("")
        pd.DataFrame(columns=["integrase_id", "passed_confident"]).to_csv(
            out_audit, sep="\t", index=False
        )
        return 0

    passed_integrase_ids: set[str] = set()

    for _, trna_row in trna_df.iterrows():
        integrase_id = str(trna_row["integrase_id"])
        contig = str(trna_row["contig"]).strip().rstrip(",")
        trna_start = int(trna_row["trna_start"])
        trna_end = int(trna_row["trna_end"])
        trna_strand = str(trna_row["trna_strand"]).strip()
        trna_len = trna_end - trna_start + 1

        if integrase_id not in int_coords:
            audit_rows.append({
                "sample": sample,
                "integrase_id": integrase_id,
                "passed_confident": False,
                "reject_reason": "integrase_not_in_hits",
            })
            continue

        int_start, int_end = int_coords[integrase_id]
        qseqid = build_qseqid(integrase_id, contig, trna_start, trna_end, trna_strand)

        raw_hits = None
        if raw_by_q is not None:
            try:
                raw_hits = raw_by_q.get_group(qseqid)
            except KeyError:
                raw_hits = pd.DataFrame()

        ie_id = next((k for k in ie_lengths if k.split(":")[0] == integrase_id), None)
        contig_seq = contig_seqs.get(contig, "")
        if not contig_seq:
            for key, val in contig_seqs.items():
                if contig in key or key in contig:
                    contig_seq = val
                    break
        attl_seq = trna_seq = ""
        site = sites.get(integrase_id)
        if site is not None and ie_id in ie_seqs:
            seq = ie_seqs[ie_id]
            attl_seq = seq[int(site["attL_start"]): int(site["attL_end"])]
            trna_seq = seq[int(site["attR_start"]): int(site["attR_end"])]
        audit = evaluate_ie_candidate(
            integrase_id=integrase_id,
            contig=contig,
            trna_start=trna_start,
            trna_end=trna_end,
            trna_strand=trna_strand,
            trna_len=trna_len,
            integrase_start=int_start,
            integrase_end=int_end,
            ie_id=ie_id,
            ie_len_nt=ie_lengths.get(ie_id) if ie_id else None,
            ie_has_n=ie_has_n.get(ie_id, True) if ie_id else True,
            raw_hits=raw_hits,
            attl_seq=attl_seq,
            trna_seq=trna_seq,
            contig_seq=contig_seq,
            cds_by_contig=cds_by_contig,
            thresholds=thresholds,
        )
        audit["sample"] = sample
        audit["qseqid"] = qseqid
        audit_rows.append(audit)
        if audit["passed_confident"]:
            passed_integrase_ids.add(integrase_id)

    audit_df = pd.DataFrame(audit_rows)
    out_audit.parent.mkdir(parents=True, exist_ok=True)
    audit_df.to_csv(out_audit, sep="\t", index=False)

    passed_records = []
    if mge_fa_path.is_file() and mge_fa_path.stat().st_size > 0:
        for rec in SeqIO.parse(mge_fa_path, "fasta"):
            if rec.id.split(":")[0] in passed_integrase_ids:
                passed_records.append(rec)
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    with out_fa.open("w") as fh:
        SeqIO.write(passed_records, fh, "fasta")

    passed_gbk = []
    if mge_gbk_path.is_file() and mge_gbk_path.stat().st_size > 0:
        for rec in SeqIO.parse(mge_gbk_path, "genbank"):
            locus_key = rec.id.split(":")[0] if ":" in rec.id else rec.id
            if locus_key in passed_integrase_ids:
                passed_gbk.append(rec)
    with out_gbk.open("w") as fh:
        SeqIO.write(passed_gbk, fh, "genbank")

    n_passed = len(passed_integrase_ids)
    logger.info(f"{sample}: {n_passed}/{len(trna_df)} confident IE(s) -> {out_fa.name}")
    return n_passed


def main() -> None:
    """Command-line entry point for per-sample confident IE filtering."""
    parser = argparse.ArgumentParser(
        description="Filter integrative element candidates to a confident set."
    )
    parser.add_argument("--sample", required=True, help="Sample identifier.")
    parser.add_argument("--trna", required=True, help="Path to integrase_trna.tsv.")
    parser.add_argument(
        "--integrase-hits",
        required=True,
        help="Path to integrase_hits_summary.tsv.",
    )
    parser.add_argument("--blast-raw", required=True, help="Path to mge_blast_raw.tsv.")
    parser.add_argument("--mge-fa", required=True, help="Path to mge_region.fa.")
    parser.add_argument("--mge-gbk", required=True, help="Path to mge_annotated.gbk.")
    parser.add_argument("--orfs-gff", required=True, help="Path to orfs.gff.")
    parser.add_argument("--fasta", required=True, help="Path to the assembly FASTA.")
    parser.add_argument("--config", default="ie_finder_config.yaml", help="Pipeline config.")
    parser.add_argument("--out-fa", required=True, help="Output ie_confident.fa path.")
    parser.add_argument("--out-gbk", required=True, help="Output ie_confident.gbk path.")
    parser.add_argument("--out-audit", required=True, help="Output ie_filter_audit.tsv path.")
    args = parser.parse_args()

    thresholds = load_thresholds(Path(args.config))
    filter_sample(
        sample=args.sample,
        trna_path=Path(args.trna),
        integrase_hits_path=Path(args.integrase_hits),
        raw_blast_path=Path(args.blast_raw),
        mge_fa_path=Path(args.mge_fa),
        mge_gbk_path=Path(args.mge_gbk),
        orfs_gff_path=Path(args.orfs_gff),
        fasta_path=Path(args.fasta),
        out_fa=Path(args.out_fa),
        out_gbk=Path(args.out_gbk),
        out_audit=Path(args.out_audit),
        thresholds=thresholds,
    )


if __name__ == "__main__":
    main()
