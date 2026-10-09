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

from extract_trna_region import extract_trna_sequence
from genome import contig_record, read_genome
from v3ps_filters import FilterThresholds, evaluate_ie_candidate

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


class ContigLookupTest(unittest.TestCase):
    def test_lookup_is_exact(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            path = tmp / "genome.fna"
            SeqIO.write(
                [SeqRecord(Seq("ACGT"), id="11", description=""),
                 SeqRecord(Seq("TTTT"), id="1", description="")],
                path,
                "fasta",
            )
            records = read_genome(path)
        self.assertEqual(str(contig_record(records, "1").seq), "TTTT")
        with self.assertRaises(KeyError):
            contig_record(records, "2")


TRNA = "GGGCCCTTTA"


def _evaluate(contig_seq: str) -> dict:
    # tRNA at 101-110 (+), its copy at 171-180, BLAST window starting at the tRNA end.
    hits = pd.DataFrame([{
        "qseqid": "1_5:1:101-110:+", "sseqid": "subject", "pident": 100.0,
        "length": 10, "mismatch": 0, "gapopen": 0, "qstart": 1, "qend": 10,
        "sstart": 62, "send": 71, "evalue": 1e-3, "bitscore": 20.0, "wstart": 110,
    }])
    return evaluate_ie_candidate(
        integrase_id="1_5", contig="1", trna_start=101, trna_end=110, trna_strand="+",
        trna_len=10, integrase_start=120, integrase_end=165, raw_hits=hits,
        contig_seq=contig_seq, cds_by_contig={}, thresholds=FilterThresholds(),
    )


class ElementSpanTest(unittest.TestCase):
    def test_span_runs_from_attr_to_attl_on_the_assembly(self):
        row = _evaluate("C" * 100 + TRNA + "C" * 60 + TRNA + "C" * 50)
        self.assertEqual(row["ie_id"], "1_5:1:101-180")
        self.assertEqual(row["ie_len_nt"], 80)
        self.assertEqual(row["exact_anchored_bp"], 10)
        self.assertFalse(row["ie_has_ambiguous_n"])

    def test_n_inside_the_element_is_flagged(self):
        row = _evaluate("C" * 100 + TRNA + "C" * 30 + "N" + "C" * 29 + TRNA + "C" * 50)
        self.assertTrue(row["ie_has_ambiguous_n"])


if __name__ == "__main__":
    unittest.main()
