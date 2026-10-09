#!/usr/bin/env python3
"""Write attL/attR on the original assembly, as GFF3 and GenBank.

This script reads the audit table, where attL and the tRNA are in 1-based
inclusive contig coordinates, and writes features on the assembly itself.

Default output is the confident set (the manuscript cascade). ``--all-candidates``
also writes candidates that have an attL coordinate but failed a later filter.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger("export_genome_features")

SOURCE = "IE_finder"
TSV_COLUMNS = [
    "sample",
    "integrase_id",
    "contig",
    "ie_start",
    "ie_end",
    "attL_start",
    "attL_end",
    "attL_strand",
    "attR_start",
    "attR_end",
    "attR_strand",
    "attR_type",
    "integrase_start",
    "integrase_end",
    "integrase_strand",
    "confidence",
]


def as_bool(value: object) -> bool:
    """Interpret audit flags written by pandas or as text."""
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes"}


def gff_strand(value: object) -> str:
    """Map Prodigal and audit strand labels onto a GFF3 strand column."""
    text = str(value).strip()
    if text in {"+", "1", "+1"}:
        return "+"
    if text in {"-", "-1"}:
        return "-"
    return "."


def strand_int(value: object) -> int | None:
    """Map a strand label onto a Biopython strand integer."""
    text = gff_strand(value)
    if text == "+":
        return 1
    if text == "-":
        return -1
    return None


def gff_escape(text: object) -> str:
    """Escape a GFF3 attribute value."""
    return (
        str(text)
        .replace("%", "%25")
        .replace(";", "%3B")
        .replace("=", "%3D")
        .replace(",", "%2C")
        .replace("\t", "%09")
        .replace("\n", "%0A")
        .replace("\r", "%0D")
    )


def ordered_span(start: int, end: int) -> tuple[int, int]:
    """Return a 1-based inclusive interval with the lower coordinate first."""
    return (start, end) if start <= end else (end, start)


def _read_table(path: Path) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", dtype={"contig": str})
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _int_field(row: pd.Series, column: str) -> int | None:
    if column not in row.index or pd.isna(row[column]):
        return None
    try:
        value = int(row[column])
    except (TypeError, ValueError):
        return None
    return value if value > 0 else None


def load_integrases(path: Path) -> dict[str, pd.Series]:
    """Index integrase nucleotide coordinates by ORF id."""
    frame = _read_table(path)
    if frame.empty or "orf_id" not in frame.columns:
        return {}
    frame["orf_id"] = frame["orf_id"].astype(str)
    return {str(row["orf_id"]): row for _, row in frame.drop_duplicates("orf_id").iterrows()}


def load_trna_types(path: Path | None) -> dict[tuple[str, int, int], str]:
    """Map (integrase id, tRNA start, tRNA end) to the Aragorn amino-acid label."""
    if path is None:
        return {}
    frame = _read_table(path)
    if frame.empty or "integrase_id" not in frame.columns:
        return {}
    type_col = "tRNA_type" if "tRNA_type" in frame.columns else None
    if type_col is None:
        return {}
    out: dict[tuple[str, int, int], str] = {}
    for _, row in frame.iterrows():
        start, end = _int_field(row, "trna_start"), _int_field(row, "trna_end")
        if start is None or end is None:
            continue
        lo, hi = ordered_span(start, end)
        out[(str(row["integrase_id"]), lo, hi)] = str(row[type_col])
    return out


def select_rows(audit: pd.DataFrame, all_candidates: bool) -> pd.DataFrame:
    """Keep confident elements, plus attL-bearing candidates when requested."""
    if audit.empty or "integrase_id" not in audit.columns:
        return pd.DataFrame()
    kept: list[pd.Series] = []
    for _, row in audit.iterrows():
        passed = "passed_confident" in row.index and as_bool(row["passed_confident"])
        has_attl = _int_field(row, "attL_abs_lo") is not None and _int_field(row, "attL_abs_hi") is not None
        if passed or (all_candidates and has_attl):
            kept.append(row)
    if not kept:
        return pd.DataFrame()
    return pd.DataFrame(kept)


def build_elements(
    audit: pd.DataFrame,
    integrases: dict[str, pd.Series],
    trna_types: dict[tuple[str, int, int], str],
    records: dict[str, SeqRecord],
    *,
    sample: str,
    all_candidates: bool,
) -> list[dict]:
    """Collect one feature bundle per element, in contig coordinates."""
    chosen = select_rows(audit, all_candidates)
    elements: list[dict] = []
    for _, row in chosen.iterrows():
        contig = str(row.get("contig", "")).strip()
        if contig not in records:
            raise SystemExit(f"Contig {contig!r} from the audit is not in the FASTA")
        attl_lo, attl_hi = _int_field(row, "attL_abs_lo"), _int_field(row, "attL_abs_hi")
        trna_lo, trna_hi = _int_field(row, "trna_start"), _int_field(row, "trna_end")
        if None in (attl_lo, attl_hi, trna_lo, trna_hi):
            logger.warning(f"Skipping {row.get('integrase_id')}: incomplete att coordinates")
            continue
        attl_start, attl_end = ordered_span(attl_lo, attl_hi)
        attr_start, attr_end = ordered_span(trna_lo, trna_hi)
        integrase_id = str(row["integrase_id"])
        hit = integrases.get(integrase_id)
        int_start = int_end = None
        int_strand = "."
        if hit is not None:
            int_start_raw, int_end_raw = _int_field(hit, "start"), _int_field(hit, "end")
            if int_start_raw is not None and int_end_raw is not None:
                int_start, int_end = ordered_span(int_start_raw, int_end_raw)
            if "strand" in hit.index:
                int_strand = gff_strand(hit["strand"])
        bounds = [attl_start, attl_end, attr_start, attr_end]
        if int_start is not None and int_end is not None:
            bounds.extend([int_start, int_end])
        span = (min(bounds), max(bounds))
        passed = "passed_confident" in row.index and as_bool(row["passed_confident"])
        elements.append({
            "sample": sample,
            "integrase_id": integrase_id,
            "contig": contig,
            "ie_start": span[0],
            "ie_end": span[1],
            "attL_start": attl_start,
            "attL_end": attl_end,
            "attL_strand": gff_strand(row.get("attL_strand", ".")),
            "attR_start": attr_start,
            "attR_end": attr_end,
            "attR_strand": gff_strand(row.get("trna_strand", ".")),
            "attR_type": trna_types.get((integrase_id, attr_start, attr_end), ""),
            "integrase_start": int_start or 0,
            "integrase_end": int_end or 0,
            "integrase_strand": int_strand,
            "confidence": "confident" if passed else "candidate",
        })
    elements.sort(key=lambda item: (item["contig"], item["ie_start"], item["integrase_id"]))
    for index, item in enumerate(elements, start=1):
        item["name"] = f"IE{index}"
    return elements


def _gff_line(contig: str, feature_type: str, start: int, end: int, strand: str, attributes: str) -> str:
    return "\t".join([
        contig, SOURCE, feature_type, str(start), str(end), ".", strand, ".", attributes,
    ])


def write_gff3(path: Path, elements: list[dict], records: dict[str, SeqRecord]) -> None:
    """Write a GFF3 file whose sequence ids are the FASTA record ids."""
    lines = ["##gff-version 3"]
    seen: list[str] = []
    for item in elements:
        if item["contig"] not in seen:
            seen.append(item["contig"])
    for contig in seen:
        lines.append(f"##sequence-region {contig} 1 {len(records[contig].seq)}")
    for item in elements:
        name = item["name"]
        parent = (
            f"ID={name};Name={name};integrase_id={gff_escape(item['integrase_id'])};"
            f"confidence={item['confidence']}"
        )
        lines.append(_gff_line(
            item["contig"], "mobile_genetic_element",
            item["ie_start"], item["ie_end"], ".", parent,
        ))
        lines.append(_gff_line(
            item["contig"], "attachment_site",
            item["attL_start"], item["attL_end"], item["attL_strand"],
            f"ID={name}_attL;Parent={name};Note=attL",
        ))
        attr_attr = f"ID={name}_attR;Parent={name};Note=attR"
        if item["attR_type"]:
            attr_attr += f";product={gff_escape(item['attR_type'])}"
        lines.append(_gff_line(
            item["contig"], "attachment_site",
            item["attR_start"], item["attR_end"], item["attR_strand"], attr_attr,
        ))
        if item["integrase_start"] and item["integrase_end"]:
            lines.append(_gff_line(
                item["contig"], "CDS",
                item["integrase_start"], item["integrase_end"], item["integrase_strand"],
                f"ID={name}_integrase;Parent={name};Note=integrase;product=tyrosine integrase",
            ))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def _gb_feature(feature_type: str, start: int, end: int, strand: str, qualifiers: dict[str, list[str]]) -> SeqFeature:
    return SeqFeature(
        FeatureLocation(start - 1, end, strand=strand_int(strand)),
        type=feature_type,
        qualifiers=qualifiers,
    )


def write_genbank(path: Path, elements: list[dict], records: dict[str, SeqRecord]) -> None:
    """Write contigs that carry an element, with features in assembly coordinates."""
    by_contig: dict[str, list[dict]] = {}
    for item in elements:
        by_contig.setdefault(item["contig"], []).append(item)
    annotated: list[SeqRecord] = []
    for contig, items in by_contig.items():
        record = records[contig]
        features: list[SeqFeature] = []
        for item in items:
            name = item["name"]
            features.append(_gb_feature(
                "mobile_element", item["ie_start"], item["ie_end"], ".",
                {
                    "mobile_element_type": ["integrative element"],
                    "note": [name, item["confidence"]],
                    "label": [name],
                },
            ))
            features.append(_gb_feature(
                "misc_feature", item["attL_start"], item["attL_end"], item["attL_strand"],
                {"note": ["attL"], "label": [f"{name} attL"]},
            ))
            attr_notes = ["attR"]
            if item["attR_type"]:
                attr_notes.append(item["attR_type"])
            features.append(_gb_feature(
                "misc_feature", item["attR_start"], item["attR_end"], item["attR_strand"],
                {"note": attr_notes, "label": [f"{name} attR"]},
            ))
            if item["integrase_start"] and item["integrase_end"]:
                features.append(_gb_feature(
                    "CDS", item["integrase_start"], item["integrase_end"], item["integrase_strand"],
                    {"note": ["integrase"], "product": ["tyrosine integrase"], "label": [f"{name} integrase"]},
                ))
        record.features = features
        record.annotations["molecule_type"] = "DNA"
        annotated.append(record)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        SeqIO.write(annotated, handle, "genbank")


def write_table(path: Path, elements: list[dict]) -> None:
    """Write a coordinate table, one element per row."""
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(elements, columns=TSV_COLUMNS) if elements else pd.DataFrame(columns=TSV_COLUMNS)
    frame.to_csv(path, sep="\t", index=False)


def export_genome_features(
    *,
    fasta: Path,
    audit: Path,
    integrases: Path,
    trna: Path | None,
    out_gff3: Path,
    out_gbk: Path,
    out_tsv: Path,
    sample: str,
    all_candidates: bool,
) -> int:
    """Write GFF3, GenBank, and a coordinate table. Returns the element count."""
    records = {rec.id: rec for rec in SeqIO.parse(fasta, "fasta")}
    elements = build_elements(
        _read_table(audit),
        load_integrases(integrases),
        load_trna_types(trna),
        records,
        sample=sample,
        all_candidates=all_candidates,
    )
    write_gff3(out_gff3, elements, records)
    write_genbank(out_gbk, elements, records)
    write_table(out_tsv, elements)
    logger.info(f"{sample}: {len(elements)} element(s) -> {out_gff3.name}")
    return len(elements)


def main() -> None:
    """Parse arguments and write the three output files."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(
        description="Write attL/attR features in the coordinates of the input assembly"
    )
    parser.add_argument("--fasta", required=True, type=Path, help="Input assembly FASTA")
    parser.add_argument("--audit", required=True, type=Path, help="ie_filter_audit.tsv")
    parser.add_argument("--integrases", required=True, type=Path, help="integrase_hits_summary.tsv")
    parser.add_argument("--trna", type=Path, default=None, help="integrase_trna.tsv")
    parser.add_argument("--out-gff3", required=True, type=Path)
    parser.add_argument("--out-gbk", required=True, type=Path)
    parser.add_argument("--out-tsv", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument(
        "--all-candidates",
        action="store_true",
        help="Also write candidates that have attL coordinates but failed a later filter",
    )
    args = parser.parse_args()
    export_genome_features(
        fasta=args.fasta,
        audit=args.audit,
        integrases=args.integrases,
        trna=args.trna,
        out_gff3=args.out_gff3,
        out_gbk=args.out_gbk,
        out_tsv=args.out_tsv,
        sample=args.sample,
        all_candidates=args.all_candidates,
    )


if __name__ == "__main__":
    main()
