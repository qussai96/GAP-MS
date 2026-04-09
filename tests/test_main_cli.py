import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

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


if __name__ == "__main__":
    unittest.main()
