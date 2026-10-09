"""One report file keeps the audit and drops empty logs."""
import tempfile
import unittest
from pathlib import Path

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from write_report import build_report


class WriteReportTest(unittest.TestCase):
    def test_report_names_the_rejection_and_skips_empty_logs(self):
        with tempfile.TemporaryDirectory() as raw:
            sample = Path(raw)
            (sample / "integrase_hits_summary.tsv").write_text(
                "orf_id\tstart\tend\ncontig_1\t1\t900\ncontig_2\t1\t900\n"
            )
            (sample / "integrase_trna.tsv").write_text(
                "integrase_id\ttrna_start\ncontig_1\t10\ncontig_2\t20\n"
            )
            (sample / "ie_filter_audit.tsv").write_text(
                "integrase_id\tpassed_confident\treject_reason\tattL_abs_lo\tattL_abs_hi\t"
                "attL_strand\ttrna_start\ttrna_end\ttrna_strand\tintegrase_len_aa\t"
                "exact_anchored_bp\tbest_attl_len_bp\tprodigal_hit_class\n"
                "contig_1\tTrue\t\t10\t29\t+\t300\t370\t-\t301\t20\t20\tintergenic\n"
                "contig_2\tFalse\texact_run_too_short\t1\t8\t+\t40\t50\t+\t301\t4\t8\tintergenic\n"
            )
            (sample / "attachment_sites_genome.tsv").write_text(
                "integrase_id\tattL_start\tattL_end\ncontig_1\t10\t29\n"
            )
            (sample / "hmm_search.log").write_text("")
            (sample / "blast_mge.log").write_text("no strict hit for contig_2\n")

            text = build_report(sample, sample="strain", mge_finder="/mge", upstream_commit="abc")

        self.assertIn("code\t/mge", text)
        self.assertIn("sample\tstrain", text)
        self.assertIn("integrase_hits\t2", text)
        self.assertIn("confident\t1", text)
        self.assertIn("rejected\t1", text)
        self.assertIn("contig_2\trejected  exact_run_too_short", text)
        self.assertIn("## blast_mge.log", text)
        self.assertNotIn("hmm_search.log", text)
        self.assertIn("audit", text)
        self.assertIn("exact_run_too_short", text)


if __name__ == "__main__":
    unittest.main()
