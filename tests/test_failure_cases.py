import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from gapms.extract import extract_features
from gapms.gtf_utils import gtf_to_df_with_genes, read_scores_csv
from tests.helpers import create_mini_gapms_dataset, write_text


class TestFailureCases(unittest.TestCase):
    def test_gtf_parser_handles_weird_supplementary_attributes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gtf_path = Path(tmpdir) / "weird.gtf"
            write_text(
                gtf_path,
                """
                chr1\tsrc\ttranscript\t1\t9\t.\t+\t.\t### totally_weird_attributes ###
                chr1\tsrc\tCDS\t1\t9\t.\t+\t0\t???broken;;;attrs
                chr1\tsrc\tgene\t20\t30\t.\t-\t.\tID=geneX
                """,
            )

            df = gtf_to_df_with_genes(gtf_path)

            self.assertEqual(len(df), 3)
            self.assertIn("gene", set(df["Type"]))
            self.assertTrue(df["Prediction_score"].notna().all())

    def test_read_scores_csv_raises_for_missing_required_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores_path = Path(tmpdir) / "bad_scores.csv"
            write_text(
                scores_path,
                """
                id,score_value
                protA,0.9
                """,
            )

            with self.assertRaisesRegex(ValueError, "Could not find both 'Protein' and 'external_score'"):
                read_scores_csv(scores_path)

    def test_extract_features_raises_clear_error_for_empty_mapping(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = create_mini_gapms_dataset(Path(tmpdir))
            write_text(paths["mapping_file"], "Peptide\tProtein\tlocation\n")

            with self.assertRaisesRegex(ValueError, "No peptides could be mapped to proteins"):
                extract_features(
                    paths["prediction_gtf"],
                    paths["protein_fasta"],
                    paths["mapping_file"],
                    paths["output_dir"],
                    paths["scores_csv"],
                )

    def test_cli_failure_with_empty_mapping_is_logged_cleanly(self):
        repo_root = Path(__file__).resolve().parents[1]

        with tempfile.TemporaryDirectory() as tmpdir:
            paths = create_mini_gapms_dataset(Path(tmpdir))
            write_text(paths["mapping_file"], "Peptide\tProtein\tlocation\n")

            cmd = [
                sys.executable,
                "-m",
                "gapms.main",
                "-g",
                str(paths["prediction_gtf"]),
                "-f",
                str(paths["protein_fasta"]),
                "-p",
                str(paths["peptides_txt"]),
                "-m",
                str(paths["mapping_file"]),
                "-o",
                str(paths["output_dir"]),
            ]

            result = subprocess.run(cmd, cwd=repo_root, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)

            log_path = paths["output_dir"] / "log.txt"
            self.assertTrue(log_path.exists())
            self.assertIn("No peptides could be mapped to proteins", log_path.read_text())


if __name__ == "__main__":
    unittest.main()
