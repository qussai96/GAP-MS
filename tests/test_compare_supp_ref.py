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
            compare_dir = paths["output_dir"] / "Compare_to_Reference"
            self.assertTrue((compare_dir / "annotation_comparison_report.tsv").exists())
            self.assertTrue((compare_dir / "annotation_comparison_summary.tsv").exists())
            self.assertTrue((compare_dir / "Different" / "peptide_support_different_start.gtf").exists())
            self.assertTrue((compare_dir / "Different" / "peptide_support_different_start.faa").exists())

    def test_generate_annotation_report_auto_detects_tmap_in_compare_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = create_mini_gapms_dataset(Path(tmpdir))
            compare_dir = paths["output_dir"] / "Compare_to_Reference"
            compare_tmap = compare_dir / "gffcmp.supported_proteins.gtf.tmap"
            compare_tmap.write_text(paths["tmap_file"].read_text())

            report_df, summary_df = generate_annotation_report(
                output_dir=paths["output_dir"],
                supported_gtf=paths["prediction_gtf"],
                reference_gtf=paths["reference_gtf"],
                all_proteins_scores=paths["all_scores_tsv"],
                peptides_bed=paths["peptides_bed"],
                protein_fasta=paths["protein_fasta"],
            )

            self.assertEqual(len(report_df), 1)
            self.assertEqual(len(summary_df), 1)


if __name__ == "__main__":
    unittest.main()
