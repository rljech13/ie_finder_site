"""search_params.yaml: reading, checking, and listing in the report."""
import sys
import tempfile
import unittest
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from v3ps_filters import (
    PARAM_FIELDS,
    FilterThresholds,
    load_thresholds,
    read_params_file,
    thresholds_from_params,
    thresholds_to_params,
)
from write_report import build_report


class SearchParamsFileTest(unittest.TestCase):
    def test_repository_file_sets_every_parameter(self):
        params = read_params_file(ROOT / "search_params.yaml")
        self.assertEqual(list(params), list(PARAM_FIELDS))

    def test_repository_values(self):
        t = load_thresholds(ROOT / "search_params.yaml")
        self.assertEqual(t.attl_window_bp, 300000)
        self.assertEqual(t.attl_select_by, "bitscore")
        self.assertFalse(t.attl_reject_gapped)
        self.assertEqual(t.trna_max_distance_bp, 500)

    def test_code_defaults_match_the_file(self):
        self.assertEqual(load_thresholds(ROOT / "search_params.yaml"), FilterThresholds())

    def test_values_round_trip(self):
        t = load_thresholds(ROOT / "search_params.yaml")
        self.assertEqual(thresholds_from_params(thresholds_to_params(t)), t)

    def test_missing_file_is_an_error(self):
        with self.assertRaises(FileNotFoundError):
            load_thresholds(ROOT / "no_such_params.yaml")


class ParamCheckTest(unittest.TestCase):
    def _write(self, text):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = Path(tmp.name) / "strict.yaml"
        path.write_text(text)
        return path

    def test_omitted_keys_keep_defaults(self):
        t = load_thresholds(self._write("attl_window_bp: 100000\n"))
        self.assertEqual(t.attl_window_bp, 100000)
        self.assertEqual(t.attl_exact_min_bp, FilterThresholds().attl_exact_min_bp)

    def test_unknown_key_names_the_file(self):
        path = self._write("attl_window: 100000\n")
        with self.assertRaises(ValueError) as ctx:
            load_thresholds(path)
        self.assertIn(str(path), str(ctx.exception))
        self.assertIn("attl_window", str(ctx.exception))

    def test_bad_values_name_the_file(self):
        for text in ("attl_reject_gapped: maybe\n", "attl_window_bp: wide\n", "ie_min_nt: true\n"):
            path = self._write(text)
            with self.subTest(text=text), self.assertRaises(ValueError) as ctx:
                load_thresholds(path)
            self.assertIn(str(path), str(ctx.exception))

    def test_quoted_booleans_are_read(self):
        t = load_thresholds(self._write('attl_reject_gapped: "true"\nreject_ambiguous_n_ie: "no"\n'))
        self.assertTrue(t.attl_reject_gapped)
        self.assertFalse(t.reject_ambiguous_n_ie)

    def test_old_pipeline_config_is_still_read(self):
        old = {
            "paths": {"results_dir": "results"},
            "filters": {"attl_window_bp": 120000, "attl_reject_gapped": True},
            "reject_ambiguous_n_ie": True,
            "pfam_profiles": ["pfam/PF00589.hmm"],
        }
        t = load_thresholds(self._write(yaml.safe_dump(old)))
        self.assertEqual(t.attl_window_bp, 120000)
        self.assertTrue(t.attl_reject_gapped)
        self.assertTrue(t.reject_ambiguous_n_ie)


class ReportParamsTest(unittest.TestCase):
    def test_report_lists_the_parameters_used(self):
        with tempfile.TemporaryDirectory() as raw:
            tmp = Path(raw)
            used = tmp / "search_params.yaml"
            used.write_text(yaml.safe_dump(thresholds_to_params(FilterThresholds()), sort_keys=False))
            report = build_report(tmp, sample="S", params_path=used)
            without = build_report(tmp, sample="S")
        block = report.split("search_params\n", 1)[1].split("\n\n", 1)[0].splitlines()
        self.assertEqual([line.split("\t")[0] for line in block], list(PARAM_FIELDS))
        self.assertIn("attl_select_by\tbitscore", block)
        self.assertNotIn("search_params", without)


if __name__ == "__main__":
    unittest.main()
