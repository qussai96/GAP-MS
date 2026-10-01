import unittest
from unittest.mock import patch

import pandas as pd

from gapms.model_utils import get_high_confident_proteins, run_protein_classifier


class TestModelSelection(unittest.TestCase):
    def test_high_confidence_requires_every_distinct_junction_supported(self):
        rows = pd.DataFrame(
            [
                {
                    "Protein": "partial",
                    "protein_specific_peptides": 0,
                    "sequence_coverage": 0.2,
                    "N_terminal_peptides": 0,
                    "C_terminal_peptides": 0,
                    "gene_specific_peptides": 1,
                    "splice_peptides": 2,
                    "splice_sites": 2,
                    "supported_splice_sites": 1,
                },
                {
                    "Protein": "complete",
                    "protein_specific_peptides": 0,
                    "sequence_coverage": 0.2,
                    "N_terminal_peptides": 0,
                    "C_terminal_peptides": 0,
                    "gene_specific_peptides": 1,
                    "splice_peptides": 1,
                    "splice_sites": 2,
                    "supported_splice_sites": 2,
                },
            ]
        )

        high_confidence = get_high_confident_proteins(rows)

        self.assertEqual(high_confidence["Protein"].tolist(), ["complete"])

    @patch("gapms.model_utils.train_iterative_model")
    @patch("gapms.model_utils.apply_pretrained_model", side_effect=ValueError("bad model"))
    def test_pretrained_failure_does_not_silently_change_mode(self, _, mock_iterative):
        empty = pd.DataFrame()

        with self.assertRaisesRegex(ValueError, "bad model"):
            run_protein_classifier(empty, empty, empty)

        mock_iterative.assert_not_called()


if __name__ == "__main__":
    unittest.main()
