import tempfile
import unittest
from pathlib import Path

import pandas as pd

from gapms.peptide_utils import calculate_sequence_coverage, mapping_file_to_df
from tests.helpers import write_text


class TestPeptideUtils(unittest.TestCase):
    def test_mapping_file_to_df_drops_duplicates_and_unmapped(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            mapping_path = Path(tmpdir) / "mapping.tsv"
            write_text(
                mapping_path,
                """
                Peptide\tProtein\tlocation
                PEPTIDE1\tprotA\t1
                PEPTIDE1\tprotA\t1
                PEPTIDE2\tunmapped\t0
                PEPTIDE3\tprotB\t5
                """,
            )

            df = mapping_file_to_df(mapping_path)

            self.assertEqual(len(df), 2)
            self.assertSetEqual(set(df["Protein"]), {"protA", "protB"})
            self.assertSetEqual(set(df["peptide_length"]), {8})

    def test_calculate_sequence_coverage_merges_overlapping_intervals(self):
        peptide_df = pd.DataFrame(
            {
                "Protein": ["protA", "protA", "protA", "protB"],
                "pep_start": [1, 3, 8, 2],
                "pep_end": [4, 6, 10, 5],
                "prot_len": [10, 10, 10, 10],
            }
        )

        coverage = calculate_sequence_coverage(peptide_df)

        self.assertAlmostEqual(coverage["protA"], 0.9, places=4)
        self.assertAlmostEqual(coverage["protB"], 0.4, places=4)


if __name__ == "__main__":
    unittest.main()
