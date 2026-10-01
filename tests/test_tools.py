import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from gapms.tools import run_gffcompare, run_proteomapper


class TestTools(unittest.TestCase):
    @patch("gapms.tools.subprocess.run")
    def test_run_proteomapper_passes_paths_without_a_shell(self, mock_run):
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output with spaces"
            output_dir.mkdir()
            protein_fasta = Path(tmpdir) / "protein db.fa"
            peptides = Path(tmpdir) / "peptide list.txt"

            run_proteomapper(protein_fasta, peptides, output_dir)

            promast_call = mock_run.call_args_list[1]
            command = promast_call.args[0]
            self.assertIsInstance(command, list)
            self.assertEqual(command[-2:], [str(protein_fasta), str(peptides)])
            self.assertNotIn("shell", promast_call.kwargs)

    @patch("gapms.tools.subprocess.run")
    def test_run_gffcompare_writes_outputs_under_compare_dir(self, mock_run):
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        with tempfile.TemporaryDirectory() as tmpdir:
            compare_dir = Path(tmpdir) / "Compare_to_Reference"
            result = run_gffcompare("reference.gtf", "supported_proteins.gtf", compare_dir)

            self.assertTrue(compare_dir.exists())
            command = mock_run.call_args[0][0]
            self.assertEqual(command[0], "gffcompare")
            self.assertIn(str(compare_dir / "gffcmp"), command)
            self.assertEqual(result.returncode, 0)

    @patch("gapms.tools.subprocess.run")
    def test_run_gffcompare_moves_stray_outputs_into_compare_dir(self, mock_run):
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            compare_dir = tmp_path / "Compare_to_Reference"
            stray_tmap = tmp_path / "gffcmp.supported_proteins.gtf.tmap"
            stray_refmap = tmp_path / "gffcmp.supported_proteins.gtf.refmap"
            stray_tmap.write_text("dummy\n")
            stray_refmap.write_text("dummy\n")

            old_cwd = os.getcwd()
            os.chdir(tmp_path)
            try:
                run_gffcompare("reference.gtf", "supported_proteins.gtf", compare_dir)
            finally:
                os.chdir(old_cwd)

            self.assertFalse(stray_tmap.exists())
            self.assertFalse(stray_refmap.exists())
            self.assertTrue((compare_dir / "gffcmp.supported_proteins.gtf.tmap").exists())
            self.assertTrue((compare_dir / "gffcmp.supported_proteins.gtf.refmap").exists())


if __name__ == "__main__":
    unittest.main()
