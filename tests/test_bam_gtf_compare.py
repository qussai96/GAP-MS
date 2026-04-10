import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import matplotlib.pyplot as plt
import pandas as pd

from gapms.bam_search import compare_bam_support_to_input_gtf
from gapms.plotting import plot_parent_run_summary


class TestBamGtfCompare(unittest.TestCase):
    def test_plot_parent_run_summary_writes_combined_parent_png(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            prediction_dir = tmp_path / "prediction_search"
            bam_dir = tmp_path / "bam_search"
            comparisons_dir = tmp_path / "comparisons"

            for directory in [
                prediction_dir / "Figures",
                bam_dir / "Figures",
                prediction_dir / "Compare_to_Reference" / "Novel",
                bam_dir / "Compare_to_Reference" / "Novel",
                comparisons_dir,
            ]:
                directory.mkdir(parents=True, exist_ok=True)

            for png_path in [
                prediction_dir / "Figures" / "Peptides_mapping_bars.png",
                bam_dir / "Figures" / "Peptides_mapping_bars.png",
                comparisons_dir / "bam_vs_input_gtf_summary.png",
            ]:
                plt.figure(figsize=(2, 1.5))
                plt.plot([0, 1], [0, 1])
                plt.title(png_path.stem)
                plt.savefig(png_path, dpi=80, bbox_inches='tight')
                plt.close()

            (prediction_dir / "Compare_to_Reference" / "annotation_comparison_summary.tsv").write_text(
                "Category\tCount\tProtein_IDs\n"
                "peptide_support_different_splice\t3\t{a,b,c}\n"
                "peptide_support_different_start\t1\t{x}\n"
            )
            (bam_dir / "Compare_to_Reference" / "annotation_comparison_summary.tsv").write_text(
                "Category\tCount\tProtein_IDs\n"
                "peptide_support_different_stop\t2\t{y,z}\n"
            )
            (prediction_dir / "Compare_to_Reference" / "Novel" / "new_predicted_proteins_scores.tsv").write_text(
                "Protein\nnovel1\nnovel2\n"
            )
            (bam_dir / "Compare_to_Reference" / "Novel" / "new_predicted_proteins_scores.tsv").write_text(
                "Protein\nnovelA\n"
            )

            output_path = plot_parent_run_summary(tmp_path)
            self.assertTrue(output_path.exists())
            self.assertGreater(output_path.stat().st_size, 0)

    @patch("gapms.bam_search.run_gffcompare")
    def test_compare_bam_support_to_input_gtf_writes_summary_and_no_overlap_stats(self, mock_run_gffcompare):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            compare_dir = tmp_path / "compare_gtf_bam"
            compare_dir.mkdir(parents=True, exist_ok=True)

            input_gtf = tmp_path / "input.gtf"
            bam_supported_gtf = tmp_path / "bam_supported.gtf"
            bam_protein_fasta = tmp_path / "bam_supported.fa"
            tmap_file = compare_dir / "gffcmp.bam_supported.gtf.tmap"

            input_gtf.write_text(
                "chr1\tsrc\tCDS\t1\t90\t.\t+\t0\tID=cds1;Parent=pred1;protein_id=pred1;gene=gene1\n"
            )
            bam_supported_gtf.write_text(
                "chr1\ttransdecoder\tmRNA\t1\t90\t.\t+\t.\tID=pred1;Parent=GENE.pred1\n"
                "chr1\ttransdecoder\tCDS\t1\t90\t.\t+\t0\tID=cds.pred1;Parent=pred1\n"
                "chr1\ttransdecoder\tmRNA\t200\t290\t.\t+\t.\tID=bamNovel.p1;Parent=GENE.bamNovel\n"
                "chr1\ttransdecoder\tCDS\t200\t290\t.\t+\t0\tID=cds.bamNovel;Parent=bamNovel.p1\n"
            )
            bam_protein_fasta.write_text(
                ">pred1\nMPEPTIDE\n>bamNovel.p1\nMNOVELSEQ\n"
            )
            pd.DataFrame(
                [
                    {"qry_id": "pred1", "ref_id": "pred1", "class_code": "="},
                    {"qry_id": "bamNovel.p1", "ref_id": "-", "class_code": "u"},
                ]
            ).to_csv(tmap_file, sep="\t", index=False)

            summary_df, no_overlap_df = compare_bam_support_to_input_gtf(
                input_gtf=input_gtf,
                bam_supported_gtf=bam_supported_gtf,
                output_dir=compare_dir,
                bam_protein_fasta=bam_protein_fasta,
            )

            mock_run_gffcompare.assert_called_once_with(input_gtf, bam_supported_gtf, compare_dir)
            self.assertTrue((compare_dir / "bam_vs_input_gtf_summary.tsv").exists())
            self.assertTrue((compare_dir / "bam_vs_input_gtf_summary.png").exists())
            self.assertTrue((compare_dir / "bam_supported_no_gene_overlap.tsv").exists())
            self.assertTrue((compare_dir / "bam_supported_no_gene_overlap.gtf").exists())
            self.assertTrue((compare_dir / "bam_supported_no_gene_overlap.faa").exists())
            self.assertEqual(int(summary_df.loc[summary_df["class_code"] == "u", "count"].iloc[0]), 1)
            self.assertEqual(no_overlap_df["qry_id"].tolist(), ["bamNovel.p1"])


if __name__ == "__main__":
    unittest.main()
