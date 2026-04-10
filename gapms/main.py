#!/usr/bin/env python3
import warnings
import time
import traceback
warnings.filterwarnings("ignore")

import sys
import argparse
from pathlib import Path
from datetime import datetime

from gapms.extract import extract_features
from gapms.filter import filter_predictions
from gapms.tools import run_proteomapper, run_psauron, run_gffread, run_gffcompare
from gapms.peptides_to_genome import map_peptides_to_genome
from gapms.parse_gffcompare_tmap import parse_gffcompare_tmap
from gapms.compare_supp_ref import generate_annotation_report
from gapms.bam_search import (
    prepare_bam_search_inputs,
    compare_bam_support_to_input_gtf,
    report_high_potential_new_gene_candidates,
)
from gapms.plotting import plot_parent_run_summary

def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(
        description="Filter GTF predictions using mass-spec data and mapping information."
    )
    parser.add_argument("-f", "--proteins", type=Path, help="Path to the protein FASTA file (optional if --gtf and --assembly are provided)")
    parser.add_argument("-g", "--gtf", type=Path, help="Path to the prediction GTF file (required unless --bam is used)")
    parser.add_argument("-b", "--bam", type=Path, help="Optional RNA-seq BAM file. When provided, GAP-MS assembles transcripts with StringTie and predicts ORFs with TransDecoder before peptide searching.")
    parser.add_argument("-a", "--assembly", type=Path, help="Path to the genome assembly file (required if --proteins is not provided or when --bam is used)")
    parser.add_argument("-p", "--peptides", type=Path, required=True, help="Path to the peptides TXT file (required)")
    parser.add_argument("-m", "--mapping", type=Path, help="Optional precomputed peptide-to-protein mapping file")
    parser.add_argument("-s", "--scores", type=Path, help="Optional precomputed external scores CSV file with columns: 'Protein', 'external_score'")
    parser.add_argument("-c", "--compute_psauron", action="store_true", help="Compute PSAURON scores")
    parser.add_argument("-rg", "--reference_gtf", type=Path, help="Optional reference gtf file")
    parser.add_argument("-rf", "--reference_fasta", type=Path, help="Optional reference fasta file")
    parser.add_argument("-o", "--output", type=Path, help="Optional output directory")
    parser.add_argument("-i", "--iterative", action="store_true", help="Train an iterative XGBoost model from the high/low-confidence sets instead of using the pre-trained classifier")

    args = parser.parse_args()
    

    # Validate required arguments
    if not args.gtf and not args.bam:
        parser.error("Either --gtf or --bam must be specified.")

    if args.bam and not args.assembly:
        parser.error("If --bam is provided, --assembly must be specified.")

    if not args.proteins and not args.assembly:
        parser.error("If --proteins is not provided, --assembly must be specified.")

    log_file = None
    original_stdout = sys.stdout
    original_stderr = sys.stderr

    try:
        if args.output:
            output_dir = args.output
            output_dir.mkdir(parents=True, exist_ok=True)
        else:
            base_dir = args.gtf.parent.parent if args.gtf else args.bam.parent
            output_dir = base_dir / "GAPMS_Output"
            output_dir.mkdir(parents=True, exist_ok=True)
            
        prediction_dir = output_dir / "prediction_search"
        bam_dir = output_dir / "bam_search"
        comparisons_dir = output_dir / "comparisons"
        prediction_output_dir = prediction_dir if args.bam else output_dir

        # Prepare logging
        log_file_path = output_dir / "log.txt"
        log_file = open(log_file_path, "w", buffering=1)
        sys.stdout = log_file
        sys.stderr = log_file

        # Print to terminal
        print(f"\nStarting GAP-MS. Check the log file in {log_file_path}", file=sys.__stdout__)

        def _ensure_branch_dirs(branch_output_dir):
            branch_output_dir.mkdir(parents=True, exist_ok=True)
            (branch_output_dir / "Figures").mkdir(parents=True, exist_ok=True)
            (branch_output_dir / "Txt").mkdir(parents=True, exist_ok=True)
            (branch_output_dir / "Compare_to_Reference").mkdir(parents=True, exist_ok=True)
            (branch_output_dir / "Compare_to_Reference" / "Novel").mkdir(parents=True, exist_ok=True)

        def _resolve_external_scores(branch_label, protein_fasta, branch_output_dir, allow_provided_scores):
            if args.scores and allow_provided_scores:
                print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Using provided external scores file for {branch_label}: {args.scores}")
                return args.scores
            if args.compute_psauron:
                print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Computing PSAURON scores for {branch_label}...")
                return run_psauron(protein_fasta, branch_output_dir)
            if args.scores and not allow_provided_scores:
                print(f"Provided external scores are used only for `prediction_search`; `{branch_label}` is scored without them.")
            return None

        def _run_reference_comparison(branch_name, branch_gtf, branch_protein_fasta, branch_output_dir):
            if not args.reference_gtf:
                return

            print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Getting New Supported Proteins for {branch_name}....")
            cleaned_reference_gtf = branch_output_dir / "reference_cleaned.gtf"
            with open(args.reference_gtf, 'r') as infile, open(cleaned_reference_gtf, 'w') as outfile:
                removed_count = 0
                for line in infile:
                    if line.startswith('#'):
                        outfile.write(line)
                    else:
                        fields = line.rstrip('\n').split('\t')
                        if len(fields) >= 9:
                            strand = fields[6]
                            attributes = fields[8]

                            if strand not in ['+', '-']:
                                removed_count += 1
                                continue

                            id_fields = ['ID=', 'transcript_id=', 'protein_id=', 'gene_id=']
                            has_valid_id = False
                            for id_field in id_fields:
                                if id_field in attributes:
                                    id_part = [attr.strip() for attr in attributes.split(';') if attr.strip().startswith(id_field)]
                                    if id_part and id_part[0].split('=')[1]:
                                        has_valid_id = True
                                        break

                            if not has_valid_id:
                                removed_count += 1
                                continue

                            outfile.write(line)
                        else:
                            removed_count += 1

            if removed_count > 0:
                print(f"Removed {removed_count} rows with invalid strand or missing ID from reference GTF")

            branch_compare_dir = branch_output_dir / "Compare_to_Reference"
            gtf_out = branch_output_dir / "supported_proteins.gtf"
            run_gffcompare(cleaned_reference_gtf, gtf_out, branch_compare_dir)
            parse_gffcompare_tmap(branch_output_dir, branch_gtf, args.reference_gtf, args.reference_fasta)

            print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Generating Annotation Comparison Report for {branch_name}....")
            generate_annotation_report(
                branch_output_dir,
                branch_gtf,
                args.reference_gtf,
                protein_fasta=branch_protein_fasta,
            )

        def run_gapms_branch(branch_name, branch_gtf, branch_protein_fasta, branch_output_dir, mapping_file=None, external_scores_csv=None):
            _ensure_branch_dirs(branch_output_dir)

            if mapping_file:
                mapping = mapping_file
                print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Using provided mapping file for {branch_name}: {mapping}")
            else:
                print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Mapping peptides to proteins with Proteomapper for {branch_name}....")
                mapping = run_proteomapper(branch_protein_fasta, args.peptides, branch_output_dir)
                if args.gtf and args.bam and Path(mapping).parent != branch_output_dir:
                    mapping = branch_output_dir / Path(mapping).name

            print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Extracting Features for {branch_name}....")
            extract_features(branch_gtf, branch_protein_fasta, mapping, branch_output_dir, external_scores_csv)

            print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Running GAPMS for {branch_name}....")
            filter_predictions(
                branch_gtf,
                branch_protein_fasta,
                branch_output_dir,
                plot_external_scores_enabled=bool(external_scores_csv),
                use_iterative_training=args.iterative,
            )

            print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Mapping peptides to genome coordinates for {branch_name}....")
            map_peptides_to_genome(branch_gtf, branch_protein_fasta, mapping, branch_output_dir)

            _run_reference_comparison(branch_name, branch_gtf, branch_protein_fasta, branch_output_dir)
            return {
                "name": branch_name,
                "gtf": branch_gtf,
                "protein_fasta": branch_protein_fasta,
                "mapping": mapping,
                "output_dir": branch_output_dir,
            }

        branch_results = {}

        if args.gtf:
            prediction_output_dir.mkdir(parents=True, exist_ok=True)
            if args.proteins:
                prediction_protein_fasta = args.proteins
            else:
                print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Generating protein FASTA from GTF and genome assembly using Gffread...")
                prediction_protein_fasta = run_gffread(args.assembly, args.gtf, prediction_output_dir, "predictions")

            prediction_external_scores_csv = _resolve_external_scores(
                "prediction_search",
                prediction_protein_fasta,
                prediction_output_dir,
                allow_provided_scores=True,
            )
            branch_results["prediction_search"] = run_gapms_branch(
                "prediction_search",
                args.gtf,
                prediction_protein_fasta,
                prediction_output_dir,
                mapping_file=args.mapping,
                external_scores_csv=prediction_external_scores_csv,
            )

        if args.bam:
            print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Building BAM-derived transcript and ORF search branch...")
            bam_outputs = prepare_bam_search_inputs(args.bam, args.assembly, bam_dir)
            bam_external_scores_csv = _resolve_external_scores(
                "bam_search",
                bam_outputs["protein_fasta"],
                bam_dir,
                allow_provided_scores=not args.gtf,
            )
            branch_results["bam_search"] = run_gapms_branch(
                "bam_search",
                bam_outputs["gtf"],
                bam_outputs["protein_fasta"],
                bam_dir,
                external_scores_csv=bam_external_scores_csv,
            )

        if "prediction_search" in branch_results and "bam_search" in branch_results:
            comparisons_dir.mkdir(parents=True, exist_ok=True)
            print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Comparing bam_search against prediction_search....")
            compare_bam_support_to_input_gtf(
                input_gtf=args.gtf,
                bam_supported_gtf=bam_dir / "supported_proteins.gtf",
                output_dir=comparisons_dir,
                bam_protein_fasta=branch_results["bam_search"]["protein_fasta"],
                prediction_supported_gtf=prediction_dir / "supported_proteins.gtf",
                prediction_novel_gtf=prediction_dir / "Compare_to_Reference" / "Novel" / "new_predicted_proteins.gtf",
                bam_novel_gtf=bam_dir / "Compare_to_Reference" / "Novel" / "new_predicted_proteins.gtf",
            )
            if args.reference_gtf:
                report_high_potential_new_gene_candidates(
                    prediction_novel_gtf=prediction_dir / "Compare_to_Reference" / "Novel" / "new_predicted_proteins.gtf",
                    bam_novel_gtf=bam_dir / "supported_proteins.gtf",
                    output_dir=output_dir,
                    bam_tmap_file=comparisons_dir / "gffcmp.supported_proteins.gtf.tmap",
                )
            summary_figure = plot_parent_run_summary(output_dir)
            print(f"Combined parent summary figure written to: {summary_figure}")

        total_time = (time.time() - start_time) / 60
        print("\n✅ Pipeline completed", file=sys.__stdout__)
        print(f"\nTotal time run: {total_time:.2f} minutes", file=sys.__stdout__)

    except Exception as e:
        print(f"\n❌ An error occurred: {str(e)}", file=sys.__stdout__)
        print(f"Please check the log file for more details: {log_file_path}", file=sys.__stdout__)
        traceback.print_exc(file=log_file)        
        raise

    finally:
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        try:
            if log_file is not None:
                log_file.close()
        except Exception:
            pass

if __name__ == "__main__":
    main()
