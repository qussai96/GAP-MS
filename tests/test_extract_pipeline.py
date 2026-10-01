import tempfile
import unittest
from pathlib import Path

import pandas as pd

from gapms.extract import engineer_features, extract_features
from tests.helpers import create_mini_gapms_dataset, write_text


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
            evidence_columns = [
                "N_terminal_peptides", "C_terminal_peptides",
                "protein_specific_peptides", "gene_specific_peptides",
                "splice_peptides", "internal_peptides",
            ]
            self.assertTrue(
                scores_df[evidence_columns].le(scores_df["mapped_peptides"], axis=0).all().all()
            )
            prot_a = scores_df.loc[scores_df["Protein"] == "protA"].iloc[0]
            self.assertEqual(int(prot_a["mapped_peptides"]), 2)
            self.assertEqual(int(prot_a["protein_specific_peptides"]), 2)
            self.assertEqual(int(prot_a["gene_specific_peptides"]), 2)
            # Adjacent CDS records are not a splice junction.
            self.assertEqual(int(prot_a["splice_sites"]), 0)
            self.assertEqual(int(prot_a["supported_splice_sites"]), 0)
            self.assertEqual(int(prot_a["splice_peptides"]), 0)
            self.assertEqual(
                int(prot_a["splice_peptides"] + prot_a["internal_peptides"]),
                int(prot_a["mapped_peptides"]),
            )
            self.assertTrue((paths["output_dir"] / "Figures" / "Peptides_mapping_bars.png").exists())

    def test_extract_features_uses_transcript_ids_for_gffread_fasta(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            paths = create_mini_gapms_dataset(tmpdir)
            write_text(
                paths["prediction_gtf"],
                """
                chr1\tsrc\ttranscript\t1\t33\t0.95\t+\t.\tgene_id "geneA"; transcript_id "txA";
                chr1\tsrc\tCDS\t1\t18\t0.95\t+\t0\tgene_id "geneA"; transcript_id "txA"; protein_id "protA";
                chr1\tsrc\tCDS\t19\t33\t0.95\t+\t0\tgene_id "geneA"; transcript_id "txA"; protein_id "protA";
                chr1\tsrc\ttranscript\t101\t130\t0.40\t+\t.\tgene_id "geneB"; transcript_id "txB";
                chr1\tsrc\tCDS\t101\t130\t0.40\t+\t0\tgene_id "geneB"; transcript_id "txB"; protein_id "protB";
                """,
            )
            write_text(paths["protein_fasta"], ">txA\nMACDEFGHAAA\n>txB\nMNPQRSTVWYA")
            write_text(
                paths["mapping_file"],
                """
                Peptide\tProtein\tlocation
                MACD\ttxA\t1
                RSTV\ttxB\t3
                """,
            )

            extract_features(
                paths["prediction_gtf"],
                paths["protein_fasta"],
                paths["mapping_file"],
                paths["output_dir"],
            )

            scores_df = pd.read_csv(paths["output_dir"] / "all_proteins_scores.tsv", sep="\t")
            self.assertSetEqual(set(scores_df["Protein"]), {"txA", "txB"})
            self.assertTrue((scores_df["mapped_peptides"] == 1).all())

    def test_extract_features_counts_distinct_supported_junctions(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            output = root / "output"
            output.mkdir()
            gtf = root / "prediction.gtf"
            fasta = root / "proteins.faa"
            mapping = root / "mapping.tsv"
            gtf.write_text(
                'chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tgene_id "g1"; transcript_id "p1";\n'
                'chr1\tsrc\tCDS\t20\t28\t.\t+\t0\tgene_id "g1"; transcript_id "p1";\n'
                'chr1\tsrc\tCDS\t40\t48\t.\t+\t0\tgene_id "g1"; transcript_id "p1";\n'
            )
            fasta.write_text(">p1\nMACDEFGHK\n")
            mapping.write_text(
                "Peptide\tProtein\tlocation\n"
                "ACDE\tp1\t2\n"
                "CDEF\tp1\t3\n"
            )

            extract_features(gtf, fasta, mapping, output)

            scores = pd.read_csv(output / "all_proteins_scores.tsv", sep="\t")
            row = scores.iloc[0]
            self.assertEqual(int(row["splice_sites"]), 2)
            # Two peptides cross the same first junction; this is still only
            # one distinct supported junction.
            self.assertEqual(int(row["splice_peptides"]), 2)
            self.assertEqual(int(row["supported_splice_sites"]), 1)


if __name__ == "__main__":
    unittest.main()
