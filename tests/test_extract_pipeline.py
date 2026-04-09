import tempfile
import unittest
from pathlib import Path

import pandas as pd

from gapms.extract import engineer_features, extract_features
from tests.helpers import create_mini_gapms_dataset


class TestExtractPipeline(unittest.TestCase):
    def test_engineer_features_adds_expected_columns(self):
        df = pd.DataFrame(
            {
                "sequence_coverage": [0.8],
                "mapped_peptides": [4],
                "protein_length": [20],
                "protein_specific_peptides": [3],
                "gene_specific_peptides": [2],
                "splice_peptides": [1],
                "internal_peptides": [2],
                "N_terminal_peptides": [1],
                "C_terminal_peptides": [1],
            }
        )

        engineered = engineer_features(df)

        self.assertIn("mapped_peptides_per_aa", engineered.columns)
        self.assertIn("junction_evidence", engineered.columns)
        self.assertAlmostEqual(engineered.loc[0, "mapped_peptides_per_aa"], 0.2, places=6)
        self.assertEqual(int(engineered.loc[0, "has_both_term"]), 1)

    def test_extract_features_writes_scores_and_plot_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = create_mini_gapms_dataset(Path(tmpdir))

            extract_features(
                paths["prediction_gtf"],
                paths["protein_fasta"],
                paths["mapping_file"],
                paths["output_dir"],
                paths["scores_csv"],
            )

            scores_path = paths["output_dir"] / "all_proteins_scores.tsv"
            self.assertTrue(scores_path.exists())

            scores_df = pd.read_csv(scores_path, sep="\t")
            self.assertSetEqual(set(scores_df["Protein"]), {"protA", "protB"})
            prot_a_cov = float(scores_df.loc[scores_df["Protein"] == "protA", "sequence_coverage"].iloc[0])
            self.assertGreater(prot_a_cov, 0.7)
            self.assertTrue((paths["output_dir"] / "Figures" / "Peptides_mapping_bars.png").exists())


if __name__ == "__main__":
    unittest.main()
