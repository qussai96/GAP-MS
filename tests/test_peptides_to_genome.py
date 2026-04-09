import tempfile
import unittest
from pathlib import Path

from gapms.peptides_to_genome import map_peptides_to_genome
from gapms.compare_supp_ref import generate_annotation_report
from tests.helpers import create_mini_gapms_dataset, write_text


class TestPeptidesToGenomeRobustness(unittest.TestCase):
    def test_map_peptides_to_genome_handles_unmatched_peptides(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = create_mini_gapms_dataset(Path(tmpdir))
            write_text(
                paths["mapping_file"],
                """
                Peptide\tProtein\tlocation
                ZZZZ\tprotA\t1
                XXXX\tprotB\t2
                """,
            )

            map_peptides_to_genome(
                paths["prediction_gtf"],
                paths["protein_fasta"],
                paths["mapping_file"],
                paths["output_dir"],
            )

            bed_path = paths["output_dir"] / "peptides_mapped.bed"
            self.assertTrue(bed_path.exists())
            lines = bed_path.read_text().strip().splitlines()
            self.assertEqual(len(lines), 1)
            self.assertTrue(lines[0].startswith("track name=Peptides"))

    def test_generate_annotation_report_handles_missing_bed_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = create_mini_gapms_dataset(Path(tmpdir))
            missing_bed = Path(tmpdir) / "missing_peptides.bed"

            report_df, summary_df = generate_annotation_report(
                output_dir=paths["output_dir"],
                supported_gtf=paths["prediction_gtf"],
                reference_gtf=paths["reference_gtf"],
                tmap_file=paths["tmap_file"],
                all_proteins_scores=paths["all_scores_tsv"],
                peptides_bed=missing_bed,
                protein_fasta=paths["protein_fasta"],
            )

            self.assertEqual(len(report_df), 1)
            self.assertEqual(len(summary_df), 1)
            self.assertTrue((paths["output_dir"] / "Compare_to_Reference" / "annotation_comparison_report.tsv").exists())


if __name__ == "__main__":
    unittest.main()
