import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from gapms.bam_search import (
    convert_transdecoder_to_genome_gtf,
    prepare_bam_search_inputs,
    run_gffread_transcripts,
)


class TestBamSearch(unittest.TestCase):
    @patch("gapms.bam_search.subprocess.run")
    def test_prepare_bam_search_inputs_runs_expected_tools(self, mock_run):
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            bam_dir = output_dir / "bam_search"
            bam_dir.mkdir(parents=True, exist_ok=True)
            (bam_dir / "transcripts.gtf").write_text(
                "chr1\tStringTie\ttranscript\t100\t399\t1000\t+\t.\tgene_id \"STRG.1\"; transcript_id \"STRG.1.1\";\n"
                "chr1\tStringTie\texon\t100\t199\t1000\t+\t.\tgene_id \"STRG.1\"; transcript_id \"STRG.1.1\"; exon_number \"1\";\n"
                "chr1\tStringTie\texon\t300\t399\t1000\t+\t.\tgene_id \"STRG.1\"; transcript_id \"STRG.1.1\"; exon_number \"2\";\n"
            )
            (bam_dir / "transcripts.fa.transdecoder.gff3").write_text(
                "##gff-version 3\n"
                "STRG.1.1\ttransdecoder\tgene\t1\t200\t.\t+\t.\tID=GENE.STRG.1.1~~STRG.1.1.p1\n"
                "STRG.1.1\ttransdecoder\tmRNA\t1\t200\t.\t+\t.\tID=STRG.1.1.p1;Parent=GENE.STRG.1.1~~STRG.1.1.p1\n"
                "STRG.1.1\ttransdecoder\tCDS\t50\t170\t.\t+\t0\tID=cds.STRG.1.1.p1;Parent=STRG.1.1.p1\n"
            )

            result = prepare_bam_search_inputs(
                bam_path=Path("reads.bam"),
                assembly_path=Path("genome.fa"),
                output_dir=bam_dir,
            )
            self.assertEqual(result["assembled_gtf"], bam_dir / "transcripts.gtf")
            self.assertEqual(result["gtf"], bam_dir / "transcripts.fa.transdecoder.genome.gtf")
            self.assertTrue(result["gtf"].exists())
            self.assertEqual(result["transcript_fasta"], bam_dir / "transcripts.fa")
            self.assertEqual(result["protein_fasta"], bam_dir / "transcripts.fa.transdecoder.pep")

            commands = [call.args[0] for call in mock_run.call_args_list]
            self.assertEqual(commands[0][0], "stringtie")
            self.assertEqual(commands[1][0], "gffread")
            self.assertEqual(commands[2][0], "TransDecoder.LongOrfs")
            self.assertEqual(commands[3][0], "TransDecoder.Predict")

    def test_convert_transdecoder_to_genome_gtf_projects_orfs_onto_genome_coordinates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            transcripts_gtf = tmp_path / "transcripts.gtf"
            transdecoder_gff3 = tmp_path / "transcripts.fa.transdecoder.gff3"
            output_gtf = tmp_path / "transcripts.fa.transdecoder.genome.gtf"

            transcripts_gtf.write_text(
                "chr1\tStringTie\ttranscript\t100\t399\t1000\t+\t.\tgene_id \"STRG.1\"; transcript_id \"STRG.1.1\";\n"
                "chr1\tStringTie\texon\t100\t199\t1000\t+\t.\tgene_id \"STRG.1\"; transcript_id \"STRG.1.1\"; exon_number \"1\";\n"
                "chr1\tStringTie\texon\t300\t399\t1000\t+\t.\tgene_id \"STRG.1\"; transcript_id \"STRG.1.1\"; exon_number \"2\";\n"
            )
            transdecoder_gff3.write_text(
                "##gff-version 3\n"
                "STRG.1.1\ttransdecoder\tgene\t1\t200\t.\t+\t.\tID=GENE.STRG.1.1~~STRG.1.1.p1\n"
                "STRG.1.1\ttransdecoder\tmRNA\t1\t200\t.\t+\t.\tID=STRG.1.1.p1;Parent=GENE.STRG.1.1~~STRG.1.1.p1\n"
                "STRG.1.1\ttransdecoder\tCDS\t50\t170\t.\t+\t0\tID=cds.STRG.1.1.p1;Parent=STRG.1.1.p1\n"
            )

            convert_transdecoder_to_genome_gtf(transdecoder_gff3, transcripts_gtf, output_gtf)
            contents = output_gtf.read_text()

            self.assertIn("chr1\tGAPMS\ttranscript\t100\t399\t1000\t+", contents)
            self.assertIn("chr1\tGAPMS\tCDS\t149\t199", contents)
            self.assertIn("chr1\tGAPMS\tCDS\t300\t369", contents)
            self.assertIn('transcript_id "STRG.1.1.p1"', contents)
            self.assertNotIn("STRG.1.1\ttransdecoder\tCDS\t50\t170", contents)

    @patch("gapms.bam_search.subprocess.run")
    def test_run_gffread_transcripts_regenerates_stale_fai_on_coordinate_error(self, mock_run):
        mock_run.side_effect = [
            MagicMock(returncode=1, stdout="", stderr="GffObj::getSpliced() error: improper genomic coordinate 98607"),
            MagicMock(returncode=0, stdout="", stderr=""),
            MagicMock(returncode=0, stdout="", stderr=""),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            assembly_path = output_dir / "genome.fa"
            transcripts_gtf = output_dir / "transcripts.gtf"
            fai_path = Path(f"{assembly_path}.fai")

            assembly_path.write_text(">chr1\nATGC\n")
            transcripts_gtf.write_text("chr1\tsrc\ttranscript\t1\t4\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";\n")
            fai_path.write_text("stale\n")

            transcripts_fa = run_gffread_transcripts(assembly_path, transcripts_gtf, output_dir)

            self.assertEqual(transcripts_fa, output_dir / "transcripts.fa")
            self.assertFalse(fai_path.exists())
            commands = [call.args[0] for call in mock_run.call_args_list]
            self.assertEqual(commands[0][0], "gffread")
            self.assertEqual(commands[1][0], "samtools")
            self.assertEqual(commands[2][0], "gffread")


if __name__ == "__main__":
    unittest.main()
