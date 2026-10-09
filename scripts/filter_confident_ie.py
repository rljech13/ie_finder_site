#!/usr/bin/env python3
"""Apply confidence filters to integrative element candidates for one sample."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from genome import contig_record, read_genome
from logger import get_logger
from v3ps_filters import (
    FilterThresholds,
    build_qseqid,
    closest_trna_rows,
    evaluate_ie_candidate,
    load_thresholds,
    parse_orfs_gff,
)

logger = get_logger("filter_confident_ie")


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
    orfs_gff_path: Path,
    fasta_path: Path,
    out_audit: Path,
    thresholds: FilterThresholds,
) -> int:
    """Evaluate every integrase-tRNA candidate of one sample and write the audit.

    Args:
        sample: Sample identifier (genome name).
        trna_path: Path to ``integrase_trna.tsv``.
        integrase_hits_path: Path to ``integrase_hits_summary.tsv``.
        raw_blast_path: Path to ``mge_blast_raw.tsv``.
        orfs_gff_path: Path to Prodigal ``orfs.gff``.
        fasta_path: Path to the assembly FASTA. attL and the element are read from it.
        out_audit: Output path for ``ie_filter_audit.tsv``.
        thresholds: Confidence filter thresholds.

    Returns:
        Number of integrase identifiers that passed all confidence filters.
    """
    trna_df = closest_trna_rows(pd.read_csv(trna_path, sep="\t", dtype={"contig": str}))
    out_audit.parent.mkdir(parents=True, exist_ok=True)
    if trna_df.empty:
        pd.DataFrame(columns=["integrase_id", "passed_confident"]).to_csv(
            out_audit, sep="\t", index=False
        )
        return 0

    int_coords = integrase_coords(integrase_hits_path)
    raw_blast = load_raw_blast(raw_blast_path)
    raw_by_q = raw_blast.groupby("qseqid") if not raw_blast.empty else None
    cds_by_contig = parse_orfs_gff(orfs_gff_path)
    genome = read_genome(fasta_path)

    audit_rows: list[dict] = []
    n_passed = 0
    for _, trna_row in trna_df.iterrows():
        integrase_id = str(trna_row["integrase_id"])
        contig = str(trna_row["contig"]).strip()
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

        audit = evaluate_ie_candidate(
            integrase_id=integrase_id,
            contig=contig,
            trna_start=trna_start,
            trna_end=trna_end,
            trna_strand=trna_strand,
            trna_len=trna_len,
            integrase_start=int_start,
            integrase_end=int_end,
            raw_hits=raw_hits,
            contig_seq=str(contig_record(genome, contig).seq).upper(),
            cds_by_contig=cds_by_contig,
            thresholds=thresholds,
        )
        audit["sample"] = sample
        audit["qseqid"] = qseqid
        audit_rows.append(audit)
        n_passed += bool(audit["passed_confident"])

    pd.DataFrame(audit_rows).to_csv(out_audit, sep="\t", index=False)
    logger.info(f"{sample}: {n_passed}/{len(trna_df)} confident IE(s) -> {out_audit.name}")
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
    parser.add_argument("--orfs-gff", required=True, help="Path to orfs.gff.")
    parser.add_argument("--fasta", required=True, help="Path to the assembly FASTA.")
    parser.add_argument(
        "--params", "--config", dest="params", default="search_params.yaml",
        help="Search parameters (search_params.yaml).",
    )
    parser.add_argument("--out-audit", required=True, help="Output ie_filter_audit.tsv path.")
    args = parser.parse_args()

    thresholds = load_thresholds(Path(args.params))
    filter_sample(
        sample=args.sample,
        trna_path=Path(args.trna),
        integrase_hits_path=Path(args.integrase_hits),
        raw_blast_path=Path(args.blast_raw),
        orfs_gff_path=Path(args.orfs_gff),
        fasta_path=Path(args.fasta),
        out_audit=Path(args.out_audit),
        thresholds=thresholds,
    )


if __name__ == "__main__":
    main()
