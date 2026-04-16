import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import call, patch

from gapms import main as gapms_main
from tests.helpers import create_mini_gapms_dataset


class TestMainCli(unittest.TestCase):
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
            bam_dir = output_dir / "bam_search"
            generated_gtf = bam_dir / "transcripts.fa.transdecoder.gff3"
            generated_fasta = bam_dir / "transcripts.fa.transdecoder.pep"

            bam_path.write_text("bam\n")
            assembly_path.write_text(">chr1\nATGCATGCATGC\n")
            peptides_path.write_text("PEPTIDE\n")

            mock_prepare_bam_search_inputs.return_value = {
                "assembled_gtf": bam_dir / "transcripts.gtf",
                "gtf": generated_gtf,
                "transcript_fasta": bam_dir / "transcripts.fa",
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

            mock_prepare_bam_search_inputs.assert_called_once_with(bam_path, assembly_path, bam_dir)
            mock_extract_features.assert_called_once_with(generated_gtf, generated_fasta, mock_run_proteomapper.return_value, bam_dir, None)
            mock_filter_predictions.assert_called_once()
            mock_map_peptides_to_genome.assert_called_once_with(generated_gtf, generated_fasta, mock_run_proteomapper.return_value, bam_dir)

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
