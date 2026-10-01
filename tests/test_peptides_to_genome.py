import tempfile
import unittest
from pathlib import Path

import pandas as pd

from gapms.peptides_to_genome import (
    export_to_bed,
    get_protein_to_cds_mapping,
    map_codon_to_genome,
    map_peptides_to_genome,
)
from gapms.compare_supp_ref import generate_annotation_report
from tests.helpers import create_mini_gapms_dataset, write_text


class TestPeptidesToGenomeRobustness(unittest.TestCase):
    def test_bed_export_sorts_mixed_numeric_and_text_sequence_ids(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "mixed.bed"
            mappings = [
                {
                    "seqid": 1, "strand": "+", "peptide": "AAA",
                    "protein": "p1", "type": "internal", "pep_start": 1,
                    "pep_end": 3, "prot_len": 3, "previous_aa": "",
                    "genomic_regions": [(10, 18)],
                },
                {
                    "seqid": "chr2", "strand": "+", "peptide": "BBB",
                    "protein": "p2", "type": "internal", "pep_start": 1,
                    "pep_end": 3, "prot_len": 3, "previous_aa": "",
                    "genomic_regions": [(20, 28)],
                },
            ]

            export_to_bed(mappings, output)

            rows = output.read_text().splitlines()[1:]
            self.assertEqual(len(rows), 2)
            self.assertTrue(rows[0].startswith("1\t"))
            self.assertTrue(rows[1].startswith("chr2\t"))

    def test_cds_mapping_is_limited_to_peptide_mapped_proteins(self):
        gtf = pd.DataFrame(
            [
                {
                    "Type": "CDS", "Protein": "p1", "Seqid": "chr1",
                    "Start": 1, "End": 9, "Strand": "+", "Frame": "0",
                },
                {
                    "Type": "CDS", "Protein": "p2", "Seqid": "chr1",
                    "Start": 20, "End": 28, "Strand": "+", "Frame": "0",
                },
            ]
        )

        mapping = get_protein_to_cds_mapping(gtf, protein_ids={"p2"})

        self.assertEqual(set(mapping), {"p2"})

    def test_map_codon_to_genome_handles_split_codon_and_phase(self):
        cds = pd.DataFrame(
            [
                {"Start": 100, "End": 103, "Strand": "+", "Frame": "0"},
                {"Start": 200, "End": 204, "Strand": "+", "Frame": "2"},
            ]
        )

        self.assertEqual(map_codon_to_genome(2, 2, cds, "+"), [(103, 103), (200, 201)])

    def test_map_codon_to_genome_handles_reverse_strand_and_initial_phase(self):
        cds = pd.DataFrame(
            [
                {"Start": 200, "End": 204, "Strand": "-", "Frame": "2"},
                {"Start": 300, "End": 304, "Strand": "-", "Frame": "1"},
            ]
        )

        self.assertEqual(map_codon_to_genome(2, 2, cds, "-"), [(203, 204), (300, 300)])

    def test_map_codon_to_genome_does_not_extend_past_annotation(self):
        cds = pd.DataFrame(
            [{"Start": 100, "End": 108, "Strand": "+", "Frame": "1"}]
        )

        self.assertEqual(map_codon_to_genome(1, 1, cds, "+"), [(101, 103)])
        self.assertEqual(map_codon_to_genome(3, 3, cds, "+"), [])

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

    def test_mapping_keeps_same_genomic_peptide_for_each_protein(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            output = root / "output"
            gtf = root / "predictions.gtf"
            fasta = root / "proteins.faa"
            mapping = root / "mapping.tsv"
            gtf.write_text(
                'chr1\tsrc\tCDS\t1\t12\t.\t+\t0\tgene_id "g1"; transcript_id "p1";\n'
                'chr1\tsrc\tCDS\t1\t12\t.\t+\t0\tgene_id "g2"; transcript_id "p2";\n'
            )
            fasta.write_text(">p1\nMACD\n>p2\nMACD\n")
            mapping.write_text(
                "Peptide\tProtein\tlocation\n"
                "MACD\tp1\t0\n"
                "MACD\tp2\t0\n"
            )

            map_peptides_to_genome(gtf, fasta, mapping, output)

            mappings = [
                line for line in (output / "peptides_mapped.bed").read_text().splitlines()
                if not line.startswith("track")
            ]
            self.assertEqual(len(mappings), 2)
            self.assertTrue(any("protein=p1" in line for line in mappings))
            self.assertTrue(any("protein=p2" in line for line in mappings))

    def test_mapping_preserves_repeated_one_based_proteomapper_locations(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            output = root / "output"
            gtf = root / "predictions.gtf"
            fasta = root / "proteins.faa"
            mapping = root / "mapping.tsv"
            gtf.write_text(
                'chr1\tsrc\tCDS\t1\t24\t.\t+\t0\tgene_id "g1"; transcript_id "p1";\n'
            )
            fasta.write_text(">p1\nMACDMACD\n")
            mapping.write_text(
                "Peptide\tProtein\tlocation\n"
                "MACD\tp1\t1\n"
                "MACD\tp1\t5\n"
            )

            map_peptides_to_genome(gtf, fasta, mapping, output)

            mappings = [
                line for line in (output / "peptides_mapped.bed").read_text().splitlines()
                if not line.startswith("track")
            ]
            self.assertEqual(len(mappings), 2)
            self.assertTrue(any("pep_start=1" in line for line in mappings))
            self.assertTrue(any("pep_start=5" in line for line in mappings))

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

            self.assertEqual(len(report_df), 0)
            self.assertEqual(len(summary_df), 0)
            self.assertTrue((paths["output_dir"] / "Compare_to_Reference" / "annotation_comparison_report.tsv").exists())
            on_disk_report = pd.read_csv(
                paths["output_dir"] / "Compare_to_Reference" / "annotation_comparison_report.tsv",
                sep="\t",
            )
            on_disk_summary = pd.read_csv(
                paths["output_dir"] / "Compare_to_Reference" / "annotation_comparison_summary.tsv",
                sep="\t",
            )
            self.assertEqual(len(on_disk_report), 0)
            self.assertEqual(len(on_disk_summary), 0)
            self.assertIn("event_coordinates", on_disk_report.columns)
            self.assertListEqual(
                on_disk_summary.columns.tolist(),
                ["Category", "Count", "Protein_IDs"],
            )


if __name__ == "__main__":
    unittest.main()
