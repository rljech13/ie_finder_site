"""Read the assembly FASTA and look up its contigs by id.

Every contig id in the pipeline comes from the first word of a FASTA header:
Prodigal, Aragorn and Biopython all cut the header at the first space. So a
lookup is exact, and a missing id means the tables and the FASTA disagree.
"""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def read_genome(path: Path | str) -> dict[str, SeqRecord]:
    """Contig id -> record for every sequence in a FASTA file."""
    return {rec.id: rec for rec in SeqIO.parse(str(path), "fasta")}


def contig_record(records: dict[str, SeqRecord], contig: str) -> SeqRecord:
    """The record whose id is ``contig``.

    Raises:
        KeyError: when no record has that id.
    """
    try:
        return records[contig]
    except KeyError:
        raise KeyError(f"Contig {contig!r} is not in the assembly FASTA") from None
