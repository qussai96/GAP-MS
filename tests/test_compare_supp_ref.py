import tempfile
import unittest
from pathlib import Path

from gapms.compare_supp_ref import classify_difference, generate_annotation_report
from tests.helpers import create_mini_gapms_dataset


class TestCompareSuppRef(unittest.TestCase):
    def test_classify_difference_detects_start_change(self):
        pred = {"seqid": "chr1", "strand": "+", "cds_regions": [(1, 18), (19, 33)]}
        ref = {"seqid": "chr1", "strand": "+", "cds_regions": [(4, 18), (19, 33)]}

        self.assertEqual(classify_difference(pred, ref), "different_start_same_stop")

    def test_generate_annotation_report_writes_expected_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = create_mini_gapms_dataset(Path(tmpdir))

            report_df, summary_df = generate_annotation_report(
                output_dir=paths["output_dir"],
                supported_gtf=paths["prediction_gtf"],
                reference_gtf=paths["reference_gtf"],
                tmap_file=paths["tmap_file"],
                all_proteins_scores=paths["all_scores_tsv"],
                peptides_bed=paths["peptides_bed"],
                protein_fasta=paths["protein_fasta"],
            )

            self.assertEqual(len(report_df), 1)
            self.assertEqual(report_df.iloc[0]["peptide_evidence_category"], "peptide_support_different_start")
            self.assertEqual(len(summary_df), 1)
            self.assertTrue((paths["output_dir"] / "annotation_comparison_report.tsv").exists())
            self.assertTrue((paths["output_dir"] / "annotation_comparison_summary.tsv").exists())


if __name__ == "__main__":
    unittest.main()
