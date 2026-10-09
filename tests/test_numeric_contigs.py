"""Assemblies whose contigs are named with bare numbers, such as "1" and "01"."""
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from extract_mge_regions import build_region_data, extract_regions, parse_blast, parse_trna
from extract_trna_region import extract_trna_sequence

# pandas reads both names as the integer 1 unless the column is kept as text.
CONTIG_1 = "A" * 100 + "GGGCCCTTTA" + "T" * 100
CONTIG_01 = "C" * 100 + "TTTTGGGGAA" + "G" * 100


def _write_genome(tmp: Path) -> Path:
    path = tmp / "genome.fna"
    SeqIO.write(
        [
            SeqRecord(Seq(CONTIG_1), id="1", description="length=210 circular=true"),
            SeqRecord(Seq(CONTIG_01), id="01", description=""),
        ],
        path,
        "fasta",
    )
    return path


def _write_trna_table(tmp: Path, contig: str) -> Path:
    path = tmp / "integrase_trna.tsv"
    pd.DataFrame([{
        "integrase_id": "1_5",
        "model": "PF22022.2",
        "integrase_start": 120,
        "integrase_end": 200,
        "integrase_strand": -1,
        "contig": contig,
        "contig_length": 210,
        "trna_start": 101,
        "trna_end": 110,
        "trna_strand": "+",
        "tRNA_type": "tRNA-Val",
        "distance": 10,
    }]).to_csv(path, sep="\t", index=False)
    return path


class ExtractTrnaRegionTest(unittest.TestCase):
    def test_numeric_contig_name_is_extracted(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            out = tmp / "mge_query.fa"
            extract_trna_sequence(_write_genome(tmp), _write_trna_table(tmp, "1"), out)
            records = list(SeqIO.parse(out, "fasta"))
        self.assertEqual([rec.id for rec in records], ["1_5:1:101-110:+"])
        self.assertEqual(str(records[0].seq), "GGGCCCTTTA")

    def test_leading_zero_contig_is_not_read_as_another_contig(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            out = tmp / "mge_query.fa"
            extract_trna_sequence(_write_genome(tmp), _write_trna_table(tmp, "01"), out)
            records = list(SeqIO.parse(out, "fasta"))
        self.assertEqual([rec.id for rec in records], ["1_5:01:101-110:+"])
        self.assertEqual(str(records[0].seq), "TTTTGGGGAA")


class ExtractMgeRegionsTest(unittest.TestCase):
    def test_region_on_numeric_contig_is_extracted(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            genome = _write_genome(tmp)
            blast = tmp / "mge_blast.tsv"
            pd.DataFrame([{
                "integrase_id": "1_5",
                "contig": "1",
                "hit_start": 171,
                "hit_end": 180,
                "pident": 100.0,
                "length": 10,
                "evalue": 1e-3,
                "bitscore": 20.0,
            }]).to_csv(blast, sep="\t", index=False)

            region_df = build_region_data(parse_blast(blast), parse_trna(_write_trna_table(tmp, "1")))
            out = tmp / "mge_region.fa"
            extract_regions(genome, region_df, out)
            records = list(SeqIO.parse(out, "fasta"))

        self.assertEqual([rec.id for rec in records], ["1_5:1:101-180"])
        self.assertEqual(str(records[0].seq), CONTIG_1[100:180])


if __name__ == "__main__":
    unittest.main()
