"""Which 3'-anchored BLAST hit becomes the reported attL."""
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from v3ps_filters import FilterThresholds, load_thresholds, select_attl_hit, select_v3ps_strict_hit

TRNA_LEN = 75
WSTART = 1025218

# Hits of a minus-strand tRNA-Ala query from T. oshimai (integrase 1_1374).
COLUMNS = ["qstart", "qend", "sstart", "send", "length", "mismatch", "gapopen",
           "pident", "evalue", "bitscore"]
FAR_LONG_WEAK = [2, 75, 214473, 214399, 76, 13, 3, 78.947, 2.87e-08, 49.1]
NEAR_SHORT_STRONG = [30, 75, 255904, 255859, 46, 3, 0, 93.478, 2.2e-14, 69.4]
# Ignored whatever their score: not anchored at the tRNA 3' end, or wrong strand.
NOT_ANCHORED = [1, 58, 200058, 200001, 58, 0, 0, 100.0, 1e-25, 108.0]
WRONG_STRAND = [10, 75, 100000, 100065, 66, 0, 0, 100.0, 1e-30, 122.0]


def _hits(*rows):
    df = pd.DataFrame(list(rows), columns=COLUMNS)
    df["wstart"] = WSTART
    return df


class SelectAttlHitTest(unittest.TestCase):
    def setUp(self):
        self.hits = _hits(FAR_LONG_WEAK, NEAR_SHORT_STRONG, NOT_ANCHORED, WRONG_STRAND)

    def test_bitscore_takes_the_better_alignment(self):
        best = select_v3ps_strict_hit(self.hits, TRNA_LEN, "-", 3, 8, "bitscore")
        self.assertEqual((best["attL_abs_lo"], best["attL_abs_hi"]), (1281076, 1281121))
        self.assertEqual(best["best_attl_len_bp"], 46)
        self.assertAlmostEqual(best["attL_bitscore"], 69.4)
        self.assertAlmostEqual(best["attL_evalue"], 2.2e-14)
        self.assertEqual(best["n_v3ps_hits"], 2)

    def test_length_keeps_the_published_rule(self):
        best = select_v3ps_strict_hit(self.hits, TRNA_LEN, "-", 3, 8, "length")
        self.assertEqual((best["attL_abs_lo"], best["attL_abs_hi"]), (1239616, 1239690))
        self.assertEqual(best["best_attl_len_bp"], 76)

    def test_equal_bitscores_prefer_the_longer_hit(self):
        shorter = list(NEAR_SHORT_STRONG)
        longer = [28, 75, 255906, 255859, 48, 4, 0, 91.7, 2.2e-14, 69.4]
        best = select_v3ps_strict_hit(_hits(shorter, longer), TRNA_LEN, "-", 3, 8, "bitscore")
        self.assertEqual(best["best_attl_len_bp"], 48)

    def test_select_attl_hit_follows_thresholds(self):
        by_length = select_attl_hit(self.hits, TRNA_LEN, "-", FilterThresholds(attl_select_by="length"))
        by_score = select_attl_hit(self.hits, TRNA_LEN, "-", FilterThresholds())
        self.assertEqual(by_length["best_attl_len_bp"], 76)
        self.assertEqual(by_score["best_attl_len_bp"], 46)

    def test_unknown_rule_is_rejected(self):
        with self.assertRaises(ValueError):
            select_v3ps_strict_hit(self.hits, TRNA_LEN, "-", 3, 8, "evalue")


class LoadThresholdsTest(unittest.TestCase):
    def _load(self, params_yaml):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "search_params.yaml"
            path.write_text(params_yaml)
            return load_thresholds(path)

    def test_default_is_bitscore(self):
        self.assertEqual(self._load("v3ps_shift: 3\n").attl_select_by, "bitscore")

    def test_length_is_read(self):
        self.assertEqual(self._load("attl_select_by: length\n").attl_select_by, "length")

    def test_unknown_value_is_rejected(self):
        with self.assertRaises(ValueError):
            self._load("attl_select_by: longest\n")


if __name__ == "__main__":
    unittest.main()
