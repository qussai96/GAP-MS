import tempfile
import unittest
from pathlib import Path

import pandas as pd

from gapms.compare_supp_ref import (
    build_reference_event_index,
    classify_difference,
    generate_annotation_report,
    get_peptides_supporting_novel_splice_sites,
    get_terminal_supporting_peptides,
    parse_gtf_to_cds_dict,
)
from gapms.parse_gffcompare_tmap import parse_gffcompare_tmap
from tests.helpers import create_mini_gapms_dataset


class TestCompareSuppRef(unittest.TestCase):
    @staticmethod
    def _bed_location(protein, strand, blocks, pep_start, pep_end, prot_len):
        return {
            "protein": protein,
            "seqid": "chr1",
            "strand": strand,
            "blocks": blocks,
            "pep_start": pep_start,
            "pep_end": pep_end,
            "prot_len": prot_len,
        }

    def test_classify_difference_detects_start_change(self):
        pred = {"seqid": "chr1", "strand": "+", "cds_regions": [(1, 18), (19, 33)]}
        ref = {"seqid": "chr1", "strand": "+", "cds_regions": [(4, 18), (19, 33)]}

        self.assertEqual(classify_difference(pred, ref), "different_start_same_stop")

    def test_cds_parser_links_transcript_and_protein_identifiers(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gtf = Path(tmpdir) / "models.gtf"
            gtf.write_text(
                'chr1\tsrc\tCDS\t1\t30\t.\t+\t0\tgene_id "g1"; '
                'transcript_id "tx1"; protein_id "prot1";\n'
            )

            parsed = parse_gtf_to_cds_dict(gtf)

            self.assertEqual(parsed["tx1"]["protein_id"], "prot1")
            self.assertEqual(parsed["prot1"]["cds_regions"], [(1, 30)])

    def test_generate_annotation_report_writes_expected_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = create_mini_gapms_dataset(Path(tmpdir))

            report_df, summary_df = generate_annotation_report(
                output_dir=paths["output_dir"],
                supported_gtf=paths["prediction_gtf"],
                reference_gtf=paths["reference_gtf"],
                tmap_file=paths["tmap_file"],
                all_proteins_scores=paths["all_scores_tsv"],
                peptides_bed=paths["peptides_bed"],
                protein_fasta=paths["protein_fasta"],
            )

            self.assertEqual(len(report_df), 1)
            self.assertEqual(report_df.iloc[0]["peptide_evidence_category"], "peptide_support_different_start")
            self.assertEqual(len(summary_df), 1)
            compare_dir = paths["output_dir"] / "Compare_to_Reference"
            self.assertTrue((compare_dir / "annotation_comparison_report.tsv").exists())
            self.assertTrue((compare_dir / "annotation_comparison_summary.tsv").exists())
            self.assertTrue((compare_dir / "Different" / "peptide_support_different_start.gtf").exists())
            self.assertTrue((compare_dir / "Different" / "peptide_support_different_start.faa").exists())

    def test_generate_annotation_report_auto_detects_tmap_in_compare_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = create_mini_gapms_dataset(Path(tmpdir))
            compare_dir = paths["output_dir"] / "Compare_to_Reference"
            compare_tmap = compare_dir / "gffcmp.supported_proteins.gtf.tmap"
            compare_tmap.write_text(paths["tmap_file"].read_text())

            report_df, summary_df = generate_annotation_report(
                output_dir=paths["output_dir"],
                supported_gtf=paths["prediction_gtf"],
                reference_gtf=paths["reference_gtf"],
                all_proteins_scores=paths["all_scores_tsv"],
                peptides_bed=paths["peptides_bed"],
                protein_fasta=paths["protein_fasta"],
            )

            self.assertEqual(len(report_df), 1)
            self.assertEqual(len(summary_df), 1)

    def test_matching_intron_chain_does_not_hide_alternative_start(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = create_mini_gapms_dataset(Path(tmpdir))
            paths["tmap_file"].write_text(
                "qry_id\tref_id\tclass_code\n"
                "protA\trefA\t=\n"
                "protB\trefB\t=\n"
            )

            report_df, _ = generate_annotation_report(
                output_dir=paths["output_dir"],
                supported_gtf=paths["prediction_gtf"],
                reference_gtf=paths["reference_gtf"],
                tmap_file=paths["tmap_file"],
                all_proteins_scores=paths["all_scores_tsv"],
                peptides_bed=paths["peptides_bed"],
                protein_fasta=paths["protein_fasta"],
            )

            self.assertEqual(len(report_df), 1)
            self.assertEqual(
                report_df.iloc[0]["peptide_evidence_category"],
                "peptide_support_different_start",
            )

    def test_generate_annotation_report_exports_different_files_for_gtf_attributes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = create_mini_gapms_dataset(Path(tmpdir))
            paths["prediction_gtf"].write_text(
                "\n".join(
                    [
                        'chr1\tsrc\ttranscript\t1\t33\t0.95\t+\t.\tgene_id "geneA"; transcript_id "protA"; protein_id "protA";',
                        'chr1\tsrc\tCDS\t1\t18\t0.95\t+\t0\tgene_id "geneA"; transcript_id "protA"; protein_id "protA";',
                        'chr1\tsrc\tCDS\t19\t33\t0.95\t+\t0\tgene_id "geneA"; transcript_id "protA"; protein_id "protA";',
                        'chr1\tsrc\ttranscript\t101\t130\t0.40\t+\t.\tgene_id "geneB"; transcript_id "protB"; protein_id "protB";',
                        'chr1\tsrc\tCDS\t101\t130\t0.40\t+\t0\tgene_id "geneB"; transcript_id "protB"; protein_id "protB";',
                        "",
                    ]
                )
            )

            generate_annotation_report(
                output_dir=paths["output_dir"],
                supported_gtf=paths["prediction_gtf"],
                reference_gtf=paths["reference_gtf"],
                tmap_file=paths["tmap_file"],
                all_proteins_scores=paths["all_scores_tsv"],
                peptides_bed=paths["peptides_bed"],
                protein_fasta=paths["protein_fasta"],
            )

            compare_dir = paths["output_dir"] / "Compare_to_Reference"
            self.assertTrue((compare_dir / "Different" / "peptide_support_different_start.gtf").exists())
            self.assertTrue((compare_dir / "Different" / "peptide_support_different_start.faa").exists())

    def test_terminal_support_requires_exact_start_and_nontryptic_stop(self):
        pred = {"seqid": "chr1", "strand": "+", "cds_regions": [(10, 39)]}
        locations = {
            "STARTPEP": [self._bed_location("p1", "+", [(10, 33)], 1, 8, 10)],
            "SHAREDSTART": [
                self._bed_location("p1", "+", [(10, 33)], 1, 8, 10),
                self._bed_location("p2", "+", [(10, 33)], 1, 8, 10),
            ],
            "AMBIGSTART": [
                self._bed_location("p1", "+", [(10, 33)], 1, 8, 10),
                self._bed_location("p2", "+", [(100, 123)], 3, 10, 12),
            ],
            "WRONGSTART": [self._bed_location("p1", "+", [(13, 36)], 2, 9, 10)],
            "ENDPEP": [self._bed_location("p1", "+", [(16, 39)], 3, 10, 10)],
            "SHAREDEND": [
                self._bed_location("p1", "+", [(16, 39)], 3, 10, 10),
                self._bed_location("p2", "+", [(16, 39)], 3, 10, 10),
            ],
            "AMBIGEND": [
                self._bed_location("p1", "+", [(16, 39)], 3, 10, 10),
                self._bed_location("p2", "+", [(16, 39)], 3, 10, 12),
            ],
            "ENDPEPK": [self._bed_location("p1", "+", [(16, 39)], 3, 10, 10)],
        }

        self.assertEqual(
            get_terminal_supporting_peptides("p1", pred, locations, "start"),
            {"STARTPEP", "SHAREDSTART"},
        )
        self.assertEqual(
            get_terminal_supporting_peptides("p1", pred, locations, "stop"),
            {"ENDPEP", "SHAREDEND"},
        )

    def test_terminal_support_is_strand_aware(self):
        pred = {"seqid": "chr1", "strand": "-", "cds_regions": [(100, 129)]}
        locations = {
            "STARTPEP": [self._bed_location("p1", "-", [(106, 129)], 1, 8, 10)],
            "ENDPEP": [self._bed_location("p1", "-", [(100, 123)], 3, 10, 10)],
        }

        self.assertEqual(
            get_terminal_supporting_peptides("p1", pred, locations, "start"),
            {"STARTPEP"},
        )
        self.assertEqual(
            get_terminal_supporting_peptides("p1", pred, locations, "stop"),
            {"ENDPEP"},
        )

    def test_splice_support_requires_both_sites_novel_across_all_isoforms(self):
        reference = {
            "r1": {"seqid": "chr1", "strand": "+", "cds_regions": [(1, 9), (21, 30)]},
            "r2": {"seqid": "chr1", "strand": "+", "cds_regions": [(1, 12), (31, 39)]},
            "r3": {"seqid": "chr1", "strand": "+", "cds_regions": [(10, 18), (24, 30)]},
        }
        reference_index = build_reference_event_index(reference)
        peptide_locations = {
            "SPLICEPEP": [self._bed_location("p1", "+", [(7, 15), (27, 35)], 1, 6, 10)],
            "SHAREDSPLICE": [
                self._bed_location("p1", "+", [(7, 15), (27, 35)], 1, 6, 10),
                self._bed_location("p2", "+", [(7, 15), (27, 35)], 1, 6, 10),
            ],
            "AMBIGSPLICE": [
                self._bed_location("p1", "+", [(7, 15), (27, 35)], 1, 6, 10),
                self._bed_location("p2", "+", [(7, 24)], 1, 6, 10),
            ],
        }
        donor_known_locations = {
            "SPLICEPEP": [self._bed_location("p1", "+", [(7, 12), (27, 35)], 1, 6, 10)]
        }
        acceptor_known_locations = {
            "SPLICEPEP": [self._bed_location("p1", "+", [(7, 15), (24, 35)], 1, 6, 10)]
        }

        both_novel = {"seqid": "chr1", "strand": "+", "cds_regions": [(1, 15), (27, 39)]}
        donor_known = {"seqid": "chr1", "strand": "+", "cds_regions": [(1, 12), (27, 39)]}
        acceptor_known = {"seqid": "chr1", "strand": "+", "cds_regions": [(1, 15), (24, 39)]}

        self.assertEqual(
            get_peptides_supporting_novel_splice_sites(
                "p1", both_novel, reference_index, peptide_locations
            ),
            {(15, 27): {"SPLICEPEP", "SHAREDSPLICE"}},
        )
        self.assertEqual(
            get_peptides_supporting_novel_splice_sites(
                "p1", donor_known, reference_index, donor_known_locations
            ),
            {},
        )
        self.assertEqual(
            get_peptides_supporting_novel_splice_sites(
                "p1", acceptor_known, reference_index, acceptor_known_locations
            ),
            {},
        )

    def test_report_treats_three_nt_stop_difference_as_same_boundary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = create_mini_gapms_dataset(Path(tmpdir))
            # protB has a predicted stop at 130; make the reference stop 127.
            reference_text = paths["reference_gtf"].read_text().replace("101\t130", "101\t127")
            paths["reference_gtf"].write_text(reference_text)
            paths["peptides_bed"].write_text(
                "chr1\t100\t130\tENDPEP\t0\t+\t100\t130\t0\t1\t30\t0\t"
                "protein=protB\tpep_start=1\tpep_end=10\tprot_len=10\n"
            )

            report_df, _ = generate_annotation_report(
                output_dir=paths["output_dir"],
                supported_gtf=paths["prediction_gtf"],
                reference_gtf=paths["reference_gtf"],
                tmap_file=paths["tmap_file"],
                all_proteins_scores=paths["all_scores_tsv"],
                peptides_bed=paths["peptides_bed"],
                protein_fasta=paths["protein_fasta"],
            )

            stop_rows = report_df[
                report_df.get("peptide_evidence_category", "").eq("peptide_support_different_stop")
            ] if not report_df.empty else report_df
            self.assertTrue(stop_rows.empty)

    def test_report_keeps_start_and_stop_events_for_same_protein(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            output = root / "output"
            output.mkdir()
            prediction_gtf = root / "prediction.gtf"
            reference_gtf = root / "reference.gtf"
            tmap = root / "report.tmap"
            scores = root / "scores.tsv"
            bed = root / "peptides.bed"
            fasta = root / "proteins.faa"

            prediction_gtf.write_text(
                'chr1\tsrc\tCDS\t1\t36\t.\t+\t0\tgene_id "g1"; transcript_id "p1";\n'
            )
            reference_gtf.write_text(
                'chr1\tref\tCDS\t4\t30\t.\t+\t0\tgene_id "rg1"; transcript_id "r1";\n'
            )
            tmap.write_text("qry_id\tref_id\tclass_code\np1\tr1\t=\n")
            scores.write_text(
                "Protein\tmapped_peptides\tsequence_coverage\tN_terminal_peptides\tC_terminal_peptides\tsplice_peptides\n"
                "p1\t2\t1.0\t1\t1\t0\n"
            )
            bed.write_text(
                "chr1\t0\t12\tSTARTPEP\t0\t+\t0\t12\t0\t1\t12\t0\tprotein=p1\tpep_start=1\tpep_end=4\tprot_len=12\n"
                "chr1\t12\t36\tENDPEP\t0\t+\t12\t36\t0\t1\t24\t0\tprotein=p1\tpep_start=5\tpep_end=12\tprot_len=12\n"
            )
            fasta.write_text(">p1\nMAAAAAAAAAAA\n")

            report_df, _ = generate_annotation_report(
                output, prediction_gtf, reference_gtf, tmap, scores, bed, fasta
            )

            self.assertSetEqual(
                set(report_df["peptide_evidence_category"]),
                {"peptide_support_different_start", "peptide_support_different_stop"},
            )

    def test_report_includes_peptide_mapped_model_absent_from_tmap(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            output = root / "output"
            output.mkdir()
            prediction_gtf = root / "prediction.gtf"
            reference_gtf = root / "reference.gtf"
            tmap = root / "report.tmap"
            scores = root / "scores.tsv"
            bed = root / "peptides.bed"
            fasta = root / "proteins.faa"

            prediction_gtf.write_text(
                'chr1\tsrc\tCDS\t1\t30\t.\t+\t0\tgene_id "g1"; transcript_id "p1";\n'
            )
            reference_gtf.write_text(
                'chr1\tref\tCDS\t4\t30\t.\t+\t0\tgene_id "rg1"; transcript_id "r1";\n'
            )
            # p1 is intentionally absent: this models a classifier-filtered
            # query that was not included in the gffcompare input.
            tmap.write_text("qry_id\tref_id\tclass_code\n")
            scores.write_text(
                "Protein\tmapped_peptides\tsequence_coverage\n"
                "p1\t1\t0.4\n"
            )
            bed.write_text(
                "chr1\t0\t12\tMACD\t0\t+\t0\t12\t0\t1\t12\t0\t"
                "protein=p1\tpep_start=1\tpep_end=4\tprot_len=10\n"
            )
            fasta.write_text(">p1\nMACDEFGHKL\n")

            report_df, _ = generate_annotation_report(
                output, prediction_gtf, reference_gtf, tmap, scores, bed, fasta
            )

            self.assertEqual(len(report_df), 1)
            self.assertEqual(report_df.iloc[0]["predicted_id"], "p1")
            self.assertEqual(
                report_df.iloc[0]["peptide_evidence_category"],
                "peptide_support_different_start",
            )
            self.assertEqual(report_df.iloc[0]["gffcompare_class_code"], "not_evaluated")

    def test_novel_gene_requires_zero_reference_cds_overlap_for_all_isoforms(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            output = root / "output"
            compare = output / "Compare_to_Reference"
            compare.mkdir(parents=True)
            prediction_gtf = root / "prediction.gtf"
            reference_gtf = root / "reference.gtf"
            reference_fasta = root / "reference.faa"

            prediction_gtf.write_text(
                "chr1\tsrc\tCDS\t100\t129\t.\t+\t0\tID=c1;Parent=pNew;gene=geneNew\n"
                "chr1\tsrc\tCDS\t100\t129\t.\t+\t0\tID=c1b;Parent=pNew2;gene=geneNew\n"
                "chr1\tsrc\tCDS\t200\t229\t.\t+\t0\tID=c2;Parent=pMixed1;gene=geneMixed\n"
                "chr1\tsrc\tCDS\t300\t329\t.\t+\t0\tID=c3;Parent=pMixed2;gene=geneMixed\n"
                "chr1\tsrc\tCDS\t210\t225\t.\t-\t0\tID=c4;Parent=pAnti;gene=geneAnti\n"
            )
            reference_gtf.write_text(
                "chr1\tref\tCDS\t210\t239\t.\t+\t0\tgene_id \"refGene\"; transcript_id \"ref1\";\n"
            )
            reference_fasta.write_text(">ref1\nAAAAAAAAAA\n")
            # The classifier-supported FASTA intentionally omits the valid
            # novel model. Novel-gene discovery receives the complete search
            # FASTA and must therefore remain peptide-driven.
            (output / "supported_proteins.fa").write_text(
                ">pMixed1\nAAAAAAAAAA\n"
            )
            prediction_fasta = root / "all_predictions.faa"
            prediction_fasta.write_text(
                ">pNew\nAAAAAAAAAA\n>pNew2\nAAAAAAAAAA\n>pMixed1\nAAAAAAAAAA\n>pMixed2\nAAAAAAAAAA\n>pAnti\nAAAAA\n"
            )
            (output / "all_proteins_scores.tsv").write_text(
                "Protein\tmapped_peptides\n"
                "pNew\t1\n"
                "pNew2\t1\n"
                "pMixed1\t1\n"
                "pMixed2\t1\n"
                "pAnti\t1\n"
            )
            (output / "peptides_mapped.bed").write_text(
                "chr1\t99\t129\tNOVELPEP\t0\t+\t99\t129\t0\t1\t30\t0\t"
                "protein=pNew\tpep_start=1\tpep_end=10\tprot_len=10\n"
                "chr1\t99\t129\tSAMEGENE\t0\t+\t99\t129\t0\t1\t30\t0\t"
                "protein=pNew\tpep_start=1\tpep_end=10\tprot_len=10\n"
                "chr1\t99\t129\tSAMEGENE\t0\t+\t99\t129\t0\t1\t30\t0\t"
                "protein=pNew2\tpep_start=1\tpep_end=10\tprot_len=10\n"
                "chr1\t99\t129\tAMBIGPEP\t0\t+\t99\t129\t0\t1\t30\t0\t"
                "protein=pNew\tpep_start=1\tpep_end=10\tprot_len=10\n"
                "chr1\t199\t229\tAMBIGPEP\t0\t+\t199\t229\t0\t1\t30\t0\t"
                "protein=pMixed1\tpep_start=1\tpep_end=10\tprot_len=10\n"
            )
            (compare / "gffcmp.test.tmap").write_text(
                "qry_id\tref_id\tclass_code\n"
                "pNew\t-\tu\n"
                "pNew2\t-\tu\n"
                "pMixed1\tref1\to\n"
                "pMixed2\t-\tu\n"
                "pAnti\tref1\tx\n"
            )

            parse_gffcompare_tmap(
                output,
                prediction_gtf,
                reference_gtf,
                reference_fasta,
                prediction_fasta=prediction_fasta,
            )

            novel_fasta = (output / "Compare_to_Reference" / "Novel" / "new_predicted_proteins.fa").read_text()
            self.assertIn(">pNew", novel_fasta)
            self.assertIn(">pNew2", novel_fasta)
            self.assertNotIn(">pMixed1", novel_fasta)
            self.assertNotIn(">pMixed2", novel_fasta)
            self.assertNotIn(">pAnti", novel_fasta)
            novel_genes = pd.read_csv(
                output / "Compare_to_Reference" / "Novel" / "new_predicted_genes.tsv",
                sep="\t",
            )
            self.assertEqual(novel_genes["Gene"].tolist(), ["geneNew"])
            self.assertEqual(
                novel_genes.iloc[0]["Event_supporting_peptides"],
                "NOVELPEP,SAMEGENE",
            )


if __name__ == "__main__":
    unittest.main()
