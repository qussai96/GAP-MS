import tempfile
import unittest
from pathlib import Path

from gapms.gtf_utils import (
    gtf_to_df_with_genes,
    read_scores_csv,
    save_gtf_subset_all_features,
)
from tests.helpers import write_text


class TestGtfUtils(unittest.TestCase):
    def test_gtf_to_df_with_genes_parses_mixed_attribute_styles(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gtf_path = Path(tmpdir) / "mixed.gtf"
            write_text(
                gtf_path,
                """
                chr1\tsrc\ttranscript\t1\t9\t0.9\t+\t.\tgene_id \"gene1\"; transcript_id \"tx1\";
                chr1\tsrc\tCDS\t1\t9\t0.9\t+\t0\tgene_id \"gene1\"; transcript_id \"tx1\"; protein_id \"prot1\";
                chr2\tsrc\tmRNA\t20\t40\t0.7\t-\t.\tID=tx2;Parent=gene2;gene=gene2
                chr2\tsrc\tCDS\t20\t40\t0.7\t-\t0\tID=cds2;Parent=tx2;protein_id=prot2;gene=gene2
                """,
            )

            df = gtf_to_df_with_genes(gtf_path)

            self.assertIn("CDS", set(df["Type"]))
            self.assertIn("transcript", set(df["Type"]))
            self.assertIn("mRNA", set(df["Type"]))
            self.assertIn("gene1", set(df["Gene"]))
            self.assertIn("gene2", set(df["Gene"]))
            self.assertIn("prot1", set(df["Protein"]))
            self.assertIn("prot2", set(df["Protein"]))

    def test_gtf_ids_are_selected_from_the_fasta_namespace(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gtf_path = Path(tmpdir) / "ensembl.gtf"
            write_text(
                gtf_path,
                """
                chr1\tsrc\ttranscript\t1\t9\t0.9\t+\t.\tgene_id "gene1"; transcript_id "tx1";
                chr1\tsrc\tCDS\t1\t9\t0.9\t+\t0\tgene_id "gene1"; transcript_id "tx1"; protein_id "prot1";
                """,
            )

            transcript_df = gtf_to_df_with_genes(gtf_path, valid_protein_ids={"tx1"})
            protein_df = gtf_to_df_with_genes(gtf_path, valid_protein_ids={"prot1"})

            self.assertSetEqual(set(transcript_df["Protein"]), {"tx1"})
            self.assertSetEqual(set(protein_df["Protein"]), {"prot1"})

    def test_gtf_ids_match_versioned_fasta_headers(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gtf_path = Path(tmpdir) / "ensembl.gtf"
            write_text(
                gtf_path,
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tgene_id \"gene1\"; transcript_id \"tx1\"; protein_id \"prot1\";",
            )

            df = gtf_to_df_with_genes(gtf_path, valid_protein_ids={"prot1.4"})

            self.assertEqual(df.iloc[0]["Protein"], "prot1.4")

    def test_save_subset_uses_matching_transcript_namespace(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            gtf_path = tmpdir / "ensembl.gtf"
            write_text(
                gtf_path,
                """
                chr1\tsrc\ttranscript\t1\t9\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1";
                chr1\tsrc\texon\t1\t9\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1";
                chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tgene_id "gene1"; transcript_id "tx1"; protein_id "prot1";
                """,
            )

            save_gtf_subset_all_features(gtf_path, {"tx1"}, tmpdir, "subset.gtf")

            self.assertEqual(len((tmpdir / "subset.gtf").read_text().splitlines()), 3)

    def test_read_scores_csv_skips_metadata_lines(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores_path = Path(tmpdir) / "scores.csv"
            write_text(
                scores_path,
                """
                metadata line without commas
                another metadata line
                description,external_score
                protA,0.91
                protB,0.12
                """,
            )

            df = read_scores_csv(scores_path)

            self.assertEqual(list(df.columns), ["Protein", "external_score"])
            self.assertEqual(df.shape, (2, 2))
            self.assertEqual(df.iloc[0]["Protein"], "protA")
            self.assertAlmostEqual(float(df.iloc[1]["external_score"]), 0.12, places=6)


if __name__ == "__main__":
    unittest.main()
