import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import matplotlib.pyplot as plt
import pandas as pd

from gapms.bam_search import (
    _summarize_supported_overlap_from_classes,
    compare_bam_support_to_input_gtf,
    report_high_potential_new_gene_candidates,
    write_merged_supported_models,
)
from gapms.plotting import plot_parent_run_summary


class TestBamGtfCompare(unittest.TestCase):
    def test_write_merged_supported_models_writes_parent_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            prediction_gtf = tmp_path / "prediction_supported.gtf"
            bam_gtf = tmp_path / "bam_supported.gtf"
            prediction_fasta = tmp_path / "prediction_supported.fa"
            bam_fasta = tmp_path / "bam_supported.fa"

            prediction_gtf.write_text(
                "chr1\tsrc\tCDS\t1\t90\t.\t+\t0\tID=cds.predA;Parent=predA;protein_id=predA;gene=geneA\n"
            )
            bam_gtf.write_text(
                "chr1\tsrc\tCDS\t200\t290\t.\t+\t0\tgene_id \"STRG.1\"; transcript_id \"bam1\"; protein_id \"bam1\";\n"
            )
            prediction_fasta.write_text(">predA\nMPEPTIDE\n")
            bam_fasta.write_text(">bam1\nMNOVELSEQ\n>predA\nMPEPTIDE\n")

            merged_gtf, merged_fasta = write_merged_supported_models(
                prediction_supported_gtf=prediction_gtf,
                bam_supported_gtf=bam_gtf,
                prediction_supported_fasta=prediction_fasta,
                bam_supported_fasta=bam_fasta,
                output_dir=tmp_path,
            )

            self.assertTrue(merged_gtf.exists())
            self.assertTrue(merged_fasta.exists())
            merged_gtf_text = merged_gtf.read_text()
            self.assertIn("predA", merged_gtf_text)
            self.assertIn("bam1", merged_gtf_text)
            merged_fasta_text = merged_fasta.read_text()
            self.assertIn(">predA", merged_fasta_text)
            self.assertIn(">bam1", merged_fasta_text)
            self.assertEqual(merged_fasta_text.count(">predA"), 1)

    def test_supported_overlap_summary_uses_selected_gffcompare_classes(self):
        summary_df = pd.DataFrame(
            [
                {"class_code": "=", "count": 3},
                {"class_code": "j", "count": 4},
                {"class_code": "c", "count": 2},
                {"class_code": "k", "count": 1},
                {"class_code": "m", "count": 5},
                {"class_code": "n", "count": 6},
                {"class_code": "u", "count": 20},
                {"class_code": "o", "count": 7},
            ]
        )

        counts = _summarize_supported_overlap_from_classes(summary_df)

        self.assertEqual(counts["accepted_overlap"], 21)
        self.assertEqual(counts["other_classes"], 27)
        self.assertEqual(counts["total"], 48)

    def test_report_high_potential_new_gene_candidates_finds_shared_novel_loci(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            prediction_novel_gtf = tmp_path / "prediction_new.gtf"
            bam_supported_gtf = tmp_path / "bam_supported.gtf"
            bam_tmap_file = tmp_path / "gffcmp.supported_proteins.gtf.tmap"

            prediction_novel_gtf.write_text(
                "chr1\tHelixer\tCDS\t100\t180\t0.0\t+\t0\tID=predA.cds1;Parent=predA;gene=predGeneA\n"
                "chr1\tHelixer\tCDS\t220\t300\t0.0\t+\t0\tID=predA.cds2;Parent=predA;gene=predGeneA\n"
                "chr1\tHelixer\tCDS\t1000\t1100\t0.0\t+\t0\tID=predB.cds1;Parent=predB;gene=predGeneB\n"
            )
            bam_supported_gtf.write_text(
                "chr1\tGAPMS\tCDS\t800\t860\t1000\t+\t0\tgene_id \"STRG.1\"; transcript_id \"bam1\"; protein_id \"bam1\";\n"
                "chr1\tGAPMS\tCDS\t900\t980\t1000\t+\t0\tgene_id \"STRG.1\"; transcript_id \"bam1\"; protein_id \"bam1\";\n"
                "chr2\tGAPMS\tCDS\t500\t600\t1000\t+\t0\tgene_id \"STRG.2\"; transcript_id \"bam2\"; protein_id \"bam2\";\n"
            )
            pd.DataFrame(
                [
                    {"qry_id": "bam1", "ref_id": "predA", "class_code": "j"},
                    {"qry_id": "bam2", "ref_id": "-", "class_code": "u"},
                ]
            ).to_csv(bam_tmap_file, sep="\t", index=False)

            candidates_df = report_high_potential_new_gene_candidates(
                prediction_novel_gtf=prediction_novel_gtf,
                bam_novel_gtf=bam_supported_gtf,
                output_dir=tmp_path,
                bam_tmap_file=bam_tmap_file,
            )

            self.assertEqual(len(candidates_df), 1)
            self.assertEqual(candidates_df.loc[0, "prediction_protein"], "predA")
            self.assertIn("bam1", candidates_df.loc[0, "bam_proteins"])
            self.assertIn("j", candidates_df.loc[0, "bam_class_codes"])
            self.assertTrue((tmp_path / "high_potential_new_gene_candidates.tsv").exists())
            self.assertFalse((tmp_path / "high_potential_new_gene_candidates_prediction.gtf").exists())
            self.assertFalse((tmp_path / "high_potential_new_gene_candidates_bam.gtf").exists())
            self.assertFalse((tmp_path / "high_potential_new_gene_candidates_prediction.fa").exists())
            self.assertFalse((tmp_path / "high_potential_new_gene_candidates_bam.fa").exists())

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
            self.assertEqual(output_path.name, "combined_gapms_summary.png")
            self.assertTrue(output_path.exists())
            self.assertGreater(output_path.stat().st_size, 0)

    @patch("gapms.bam_search.run_gffcompare")
    def test_compare_bam_support_to_input_gtf_writes_summary_and_no_overlap_stats(self, mock_run_gffcompare):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            compare_dir = tmp_path / "compare_gtf_bam"
            compare_dir.mkdir(parents=True, exist_ok=True)

            input_gtf = tmp_path / "input.gtf"
            prediction_supported_gtf = tmp_path / "prediction_supported.gtf"
            prediction_novel_gtf = tmp_path / "prediction_novel.gtf"
            bam_supported_gtf = tmp_path / "bam_supported.gtf"
            bam_novel_gtf = tmp_path / "bam_novel.gtf"
            bam_protein_fasta = tmp_path / "bam_supported.fa"
            tmap_file = compare_dir / "gffcmp.bam_supported.gtf.tmap"

            input_gtf.write_text(
                "chr1\tsrc\tCDS\t1\t90\t.\t+\t0\tID=cds1;Parent=pred1;protein_id=pred1;gene=gene1\n"
            )
            prediction_supported_gtf.write_text(
                "chr1\tsrc\tCDS\t1\t90\t.\t+\t0\tID=cds1;Parent=pred1;protein_id=pred1;gene=gene1\n"
                "chr1\tsrc\tCDS\t400\t450\t.\t+\t0\tID=cds2;Parent=predOnly;protein_id=predOnly;gene=gene2\n"
            )
            prediction_novel_gtf.write_text(
                "chr1\tsrc\tCDS\t500\t560\t.\t+\t0\tID=cds3;Parent=predNovel;protein_id=predNovel;gene=gene3\n"
            )
            bam_supported_gtf.write_text(
                "chr1\ttransdecoder\tmRNA\t1\t90\t.\t+\t.\tID=pred1;Parent=GENE.pred1\n"
                "chr1\ttransdecoder\tCDS\t1\t90\t.\t+\t0\tID=cds.pred1;Parent=pred1\n"
                "chr1\ttransdecoder\tmRNA\t200\t290\t.\t+\t.\tID=bamNovel.p1;Parent=GENE.bamNovel\n"
                "chr1\ttransdecoder\tCDS\t200\t290\t.\t+\t0\tID=cds.bamNovel;Parent=bamNovel.p1\n"
            )
            bam_novel_gtf.write_text(
                "chr1\ttransdecoder\tmRNA\t520\t590\t.\t+\t.\tID=bamNovelShared;Parent=GENE.bamNovelShared\n"
                "chr1\ttransdecoder\tCDS\t520\t590\t.\t+\t0\tID=cds.bamNovelShared;Parent=bamNovelShared\n"
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
                prediction_supported_gtf=prediction_supported_gtf,
                prediction_novel_gtf=prediction_novel_gtf,
                bam_novel_gtf=bam_novel_gtf,
            )

            mock_run_gffcompare.assert_has_calls([
                unittest.mock.call(prediction_supported_gtf, bam_supported_gtf, compare_dir),
                unittest.mock.call(bam_supported_gtf, prediction_supported_gtf, compare_dir / "_tmp_prediction_vs_bam"),
            ])
            self.assertEqual(mock_run_gffcompare.call_count, 2)
            self.assertTrue((compare_dir / "bam_vs_prediction_supported_class_summary.tsv").exists())
            self.assertTrue((compare_dir / "bam_vs_prediction_supported_summary.png").exists())
            self.assertTrue((compare_dir / "overlapped_supported_proteins.gtf").exists())
            self.assertTrue((compare_dir / "overlapped_supported_proteins.faa").exists())
            self.assertTrue((compare_dir / "bam_supported_no_gtf_overalp.gtf").exists())
            self.assertTrue((compare_dir / "bam_supported_no_gtf_overalp.faa").exists())
            self.assertTrue((compare_dir / "gtf_supported_no_bam_overlap.gtf").exists())
            self.assertTrue((compare_dir / "gtf_supported_no_bam_overlap.faa").exists())
            self.assertTrue((compare_dir / "bam_vs_input_gtf_summary.tsv").exists())
            self.assertTrue((compare_dir / "bam_vs_input_gtf_summary.png").exists())
            self.assertTrue((compare_dir / "bam_supported_no_gene_overlap.tsv").exists())
            self.assertTrue((compare_dir / "bam_supported_no_gene_overlap.gtf").exists())
            self.assertTrue((compare_dir / "bam_supported_no_gene_overlap.faa").exists())
            self.assertEqual(int(summary_df.loc[summary_df["class_code"] == "u", "count"].iloc[0]), 1)
            self.assertEqual(no_overlap_df["qry_id"].tolist(), ["bamNovel.p1"])


if __name__ == "__main__":
    unittest.main()
