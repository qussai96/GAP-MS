import unittest

import pandas as pd

from gapms.score_filter import find_apply_score_filter, find_cutoff


class TestScoreFilter(unittest.TestCase):
    def setUp(self):
        self.df = pd.DataFrame(
            {
                "Protein": ["p1", "p2", "p3", "p4"],
                "external_score": [0.95, 0.90, 0.30, 0.05],
            }
        )
        self.high_conf = {"p1", "p2"}

    def test_find_cutoff_returns_valid_probability(self):
        cutoff = find_cutoff(self.df, self.high_conf, "external_score")
        self.assertGreaterEqual(cutoff, 0.0)
        self.assertLessEqual(cutoff, 1.0)

    def test_find_apply_score_filter_keeps_high_scores(self):
        cutoff, supported = find_apply_score_filter(self.df, self.high_conf, "external_score")
        self.assertGreaterEqual(cutoff, 0.5)
        self.assertIn("p1", supported)
        self.assertIn("p2", supported)
        self.assertNotIn("p4", supported)


if __name__ == "__main__":
    unittest.main()
