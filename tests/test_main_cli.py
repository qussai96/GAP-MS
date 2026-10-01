import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import call, patch

from gapms import main as gapms_main
from tests.helpers import create_mini_gapms_dataset


class TestMainCli(unittest.TestCase):
    def test_reference_events_use_all_peptide_mapped_models(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            output = root / "output"
            prediction_gtf = root / "predictions.gtf"
            prediction_fasta = root / "predictions.fa"
            reference_gtf = root / "reference.gtf"
            reference_fasta = root / "reference.fa"
            peptides = root / "peptides.txt"
            mapping = root / "mapping.tsv"

            prediction_gtf.write_text(
                'chr1\tsrc\tCDS\t1\t30\t.\t+\t0\tgene_id "g1"; transcript_id "P1";\n'
            )
            prediction_fasta.write_text(">P1\nMAAAAAAAAA\n>P2\nMCCCCCCCCC\n")
            reference_gtf.write_text(
                'chr1\tref\tCDS\t4\t30\t.\t+\t0\tgene_id "r1"; transcript_id "R1";\n'
            )
            reference_fasta.write_text(">R1\nAAAAAAAAA\n")
            peptides.write_text("MAAA\n")
            mapping.write_text("Peptide\tProtein\tlocation\nMAAA\tP1\t1\n")

            def fake_extract(_, __, ___, branch_output, ____):
                branch_output.mkdir(parents=True, exist_ok=True)
                (branch_output / "all_proteins_scores.tsv").write_text(
                    "Protein\tmapped_peptides\nP1\t1\nP2\t0\n"
                )

            def fake_filter(_, __, branch_output, **___):
                (branch_output / "supported_proteins.gtf").write_text("")
                (branch_output / "supported_proteins.fa").write_text("")

            def fake_save(gtf_path, proteins, branch_output, filename):
                self.assertEqual(Path(gtf_path), prediction_gtf)
                self.assertSetEqual(set(proteins), {"P1"})
                (Path(branch_output) / filename).write_text(prediction_gtf.read_text())

            argv = [
                "gapms.main",
                "-g", str(prediction_gtf),
                "-f", str(prediction_fasta),
                "-p", str(peptides),
                "-m", str(mapping),
                "-rg", str(reference_gtf),
                "-rf", str(reference_fasta),
                "-o", str(output),
            ]

            with (
                patch.object(gapms_main, "extract_features", side_effect=fake_extract),
                patch.object(gapms_main, "filter_predictions", side_effect=fake_filter),
                patch.object(gapms_main, "map_peptides_to_genome"),
                patch.object(gapms_main, "save_gtf_subset_all_features", side_effect=fake_save),
                patch.object(gapms_main, "run_gffcompare") as mock_gffcompare,
                patch.object(gapms_main, "parse_gffcompare_tmap") as mock_parse,
                patch.object(gapms_main, "generate_annotation_report") as mock_report,
                patch.object(sys, "argv", argv),
            ):
                gapms_main.main()

            query_gtf = output / "peptide_mapped_proteins_cleaned.gtf"
            mock_gffcompare.assert_called_once_with(
                output / "reference_cleaned.gtf",
                query_gtf,
                output / "Compare_to_Reference",
            )
            self.assertEqual(mock_report.call_args.args[1], prediction_gtf)
            self.assertEqual(
                mock_parse.call_args.kwargs["prediction_fasta"],
                prediction_fasta,
            )
            current_tmap = (
                output / "Compare_to_Reference"
                / "gffcmp.peptide_mapped_proteins_cleaned.gtf.tmap"
            )
            self.assertEqual(mock_parse.call_args.kwargs["tmap_file"], current_tmap)
            self.assertEqual(mock_report.call_args.kwargs["tmap_file"], current_tmap)

    def test_cli_smoke_run_with_precomputed_inputs(self):
        repo_root = Path(__file__).resolve().parents[1]

        with tempfile.TemporaryDirectory() as tmpdir:
            paths = create_mini_gapms_dataset(Path(tmpdir))

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
                "-s",
                str(paths["scores_csv"]),
                "-o",
                str(paths["output_dir"]),
            ]

            result = subprocess.run(cmd, cwd=repo_root, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, msg=result.stdout + "\n" + result.stderr)
            self.assertTrue((paths["output_dir"] / "all_proteins_scores.tsv").exists())
            self.assertTrue((paths["output_dir"] / "supported_proteins.gtf").exists())
            self.assertTrue((paths["output_dir"] / "supported_proteins.fa").exists())
            self.assertTrue((paths["output_dir"] / "peptides_mapped.bed").exists())
            self.assertTrue((paths["output_dir"] / "Txt" / "supported_proteins.txt").exists())
            log_text = (paths["output_dir"] / "log.txt").read_text()
            self.assertIn("Pipeline completed successfully", log_text)

    @patch("gapms.main.map_peptides_to_genome")
    @patch("gapms.main.filter_predictions")
    @patch("gapms.main.extract_features")
    @patch("gapms.main.run_proteomapper")
    @patch("gapms.main.prepare_bam_search_inputs")
    def test_cli_accepts_bam_input_and_uses_generated_targets(
        self,
        mock_prepare_bam_search_inputs,
        mock_run_proteomapper,
        mock_extract_features,
        mock_filter_predictions,
        mock_map_peptides_to_genome,
    ):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            bam_path = tmp_path / "reads.bam"
            assembly_path = tmp_path / "genome.fa"
            peptides_path = tmp_path / "peptides.txt"
            output_dir = tmp_path / "output"
            generated_gtf = output_dir / "transcripts.fa.transdecoder.gff3"
            generated_fasta = output_dir / "transcripts.fa.transdecoder.pep"

            bam_path.write_text("bam\n")
            assembly_path.write_text(">chr1\nATGCATGCATGC\n")
            peptides_path.write_text("PEPTIDE\n")

            mock_prepare_bam_search_inputs.return_value = {
                "assembled_gtf": output_dir / "transcripts.gtf",
                "gtf": generated_gtf,
                "transcript_fasta": output_dir / "transcripts.fa",
                "protein_fasta": generated_fasta,
            }
            mock_run_proteomapper.return_value = output_dir / "proteomapper.tsv"

            argv = [
                "gapms.main",
                "-b",
                str(bam_path),
                "-a",
                str(assembly_path),
                "-p",
                str(peptides_path),
                "-o",
                str(output_dir),
            ]

            with patch.object(sys, "argv", argv):
                gapms_main.main()

            mock_prepare_bam_search_inputs.assert_called_once_with(bam_path, assembly_path, output_dir)
            mock_extract_features.assert_called_once_with(generated_gtf, generated_fasta, mock_run_proteomapper.return_value, output_dir, None)
            mock_filter_predictions.assert_called_once()
            mock_map_peptides_to_genome.assert_called_once_with(generated_gtf, generated_fasta, mock_run_proteomapper.return_value, output_dir)

    @patch("gapms.main.compare_bam_support_to_input_gtf")
    @patch("gapms.main.write_merged_supported_models")
    @patch("gapms.main.map_peptides_to_genome")
    @patch("gapms.main.filter_predictions")
    @patch("gapms.main.extract_features")
    @patch("gapms.main.run_proteomapper")
    @patch("gapms.main.prepare_bam_search_inputs")
    @patch("gapms.main.run_gffread")
    def test_cli_creates_prediction_branch_before_gffread_when_bam_and_gtf_are_provided_without_proteins(
        self,
        mock_run_gffread,
        mock_prepare_bam_search_inputs,
        mock_run_proteomapper,
        mock_extract_features,
        mock_filter_predictions,
        mock_map_peptides_to_genome,
        mock_write_merged_supported_models,
        mock_compare_bam_support_to_input_gtf,
    ):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            bam_path = tmp_path / "reads.bam"
            assembly_path = tmp_path / "genome.fa"
            peptides_path = tmp_path / "peptides.txt"
            prediction_gtf = tmp_path / "predictions.gtf"
            output_dir = tmp_path / "output"
            prediction_dir = output_dir / "prediction_search"
            bam_dir = output_dir / "bam_search"

            bam_path.write_text("bam\n")
            assembly_path.write_text(">chr1\nATGCATGCATGC\n")
            peptides_path.write_text("PEPTIDE\n")
            prediction_gtf.write_text("chr1\tsrc\tCDS\t1\t12\t.\t+\t0\tID=cds1;Parent=prot1;protein_id=prot1;gene=gene1\n")

            def fake_run_gffread(assembly_arg, gtf_arg, branch_output_dir, label):
                self.assertEqual(branch_output_dir, prediction_dir)
                self.assertTrue(branch_output_dir.exists())
                return prediction_dir / "translated_predictions_proteins.fa"

            mock_run_gffread.side_effect = fake_run_gffread
            mock_prepare_bam_search_inputs.return_value = {
                "assembled_gtf": bam_dir / "transcripts.gtf",
                "gtf": bam_dir / "transcripts.fa.transdecoder.gff3",
                "transcript_fasta": bam_dir / "transcripts.fa",
                "protein_fasta": bam_dir / "transcripts.fa.transdecoder.pep",
            }
            mock_run_proteomapper.side_effect = [prediction_dir / "proteomapper.tsv", bam_dir / "proteomapper.tsv"]

            argv = [
                "gapms.main",
                "-g",
                str(prediction_gtf),
                "-b",
                str(bam_path),
                "-a",
                str(assembly_path),
                "-p",
                str(peptides_path),
                "-o",
                str(output_dir),
            ]

            with patch.object(sys, "argv", argv):
                gapms_main.main()

            mock_run_gffread.assert_called_once_with(assembly_path, prediction_gtf, prediction_dir, "predictions")

    @patch("gapms.main.compare_bam_support_to_input_gtf")
    @patch("gapms.main.write_merged_supported_models")
    @patch("gapms.main.map_peptides_to_genome")
    @patch("gapms.main.filter_predictions")
    @patch("gapms.main.extract_features")
    @patch("gapms.main.run_proteomapper")
    @patch("gapms.main.prepare_bam_search_inputs")
    def test_cli_runs_both_input_gtf_and_bam_branches_when_both_are_provided(
        self,
        mock_prepare_bam_search_inputs,
        mock_run_proteomapper,
        mock_extract_features,
        mock_filter_predictions,
        mock_map_peptides_to_genome,
        mock_write_merged_supported_models,
        mock_compare_bam_support_to_input_gtf,
    ):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            bam_path = tmp_path / "reads.bam"
            assembly_path = tmp_path / "genome.fa"
            peptides_path = tmp_path / "peptides.txt"
            prediction_gtf = tmp_path / "predictions.gtf"
            prediction_fasta = tmp_path / "predictions.fa"
            output_dir = tmp_path / "output"
            prediction_dir = output_dir / "prediction_search"
            bam_dir = output_dir / "bam_search"

            bam_path.write_text("bam\n")
            assembly_path.write_text(">chr1\nATGCATGCATGC\n")
            peptides_path.write_text("PEPTIDE\n")
            prediction_gtf.write_text("chr1\tsrc\tCDS\t1\t12\t.\t+\t0\tID=cds1;Parent=prot1;protein_id=prot1;gene=gene1\n")
            prediction_fasta.write_text(">prot1\nMPEPTIDE\n")

            generated_gtf = bam_dir / "transcripts.fa.transdecoder.gff3"
            generated_fasta = bam_dir / "transcripts.fa.transdecoder.pep"
            mock_prepare_bam_search_inputs.return_value = {
                "assembled_gtf": bam_dir / "transcripts.gtf",
                "gtf": generated_gtf,
                "transcript_fasta": bam_dir / "transcripts.fa",
                "protein_fasta": generated_fasta,
            }
            mock_run_proteomapper.side_effect = [output_dir / "proteomapper.tsv", bam_dir / "proteomapper.tsv"]

            argv = [
                "gapms.main",
                "-g",
                str(prediction_gtf),
                "-f",
                str(prediction_fasta),
                "-b",
                str(bam_path),
                "-a",
                str(assembly_path),
                "-p",
                str(peptides_path),
                "-o",
                str(output_dir),
            ]

            with patch.object(sys, "argv", argv):
                gapms_main.main()

            self.assertEqual(mock_run_proteomapper.call_count, 2)
            mock_run_proteomapper.assert_has_calls([
                call(prediction_fasta, peptides_path, prediction_dir),
                call(generated_fasta, peptides_path, bam_dir),
            ])
            self.assertEqual(mock_extract_features.call_count, 2)
            mock_extract_features.assert_has_calls([
                call(prediction_gtf, prediction_fasta, prediction_dir / "proteomapper.tsv", prediction_dir, None),
                call(generated_gtf, generated_fasta, bam_dir / "proteomapper.tsv", bam_dir, None),
            ])
            self.assertEqual(mock_filter_predictions.call_count, 2)
            self.assertEqual(mock_map_peptides_to_genome.call_count, 2)
            mock_compare_bam_support_to_input_gtf.assert_called_once_with(
                input_gtf=prediction_gtf,
                bam_supported_gtf=bam_dir / "supported_proteins.gtf",
                output_dir=output_dir / "comparisons",
                bam_protein_fasta=generated_fasta,
                prediction_supported_gtf=prediction_dir / "supported_proteins.gtf",
                prediction_novel_gtf=prediction_dir / "Compare_to_Reference" / "Novel" / "new_predicted_proteins.gtf",
                bam_novel_gtf=bam_dir / "Compare_to_Reference" / "Novel" / "new_predicted_proteins.gtf",
            )
            mock_write_merged_supported_models.assert_called_once_with(
                prediction_supported_gtf=prediction_dir / "supported_proteins.gtf",
                bam_supported_gtf=bam_dir / "supported_proteins.gtf",
                prediction_supported_fasta=prediction_dir / "supported_proteins.fa",
                bam_supported_fasta=bam_dir / "supported_proteins.fa",
                output_dir=output_dir,
            )


if __name__ == "__main__":
    unittest.main()
