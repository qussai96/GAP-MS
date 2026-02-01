#!/usr/bin/env python3
import warnings
import time
import traceback
warnings.filterwarnings("ignore")

import sys
import argparse
from pathlib import Path
from datetime import datetime
import subprocess

from gapms.extract import extract_features
from gapms.filter import filter_predictions
from gapms.tools import run_proteomapper, run_psauron, run_gffread, run_embeddings, run_gffcompare
from gapms.get_new_proteins import get_new_proteins
from gapms.peptides_to_genome import map_peptides_to_genome
from gapms.parse_gffcompare_tmap import parse_gffcompare_tmap

def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(
        description="Filter GTF predictions using mass-spec data and mapping information."
    )
    parser.add_argument("-f", "--proteins", type=Path, help="Path to the protein FASTA file (optional if --gtf and --assembly are provided)")
    parser.add_argument("-g", "--gtf", type=Path, required=True, help="Path to the prediction GTF file (required)")
    parser.add_argument("-a", "--assembly", type=Path, help="Path to the genome assembly file (required if --proteins is not provided)")
    parser.add_argument("-p", "--peptides", type=Path, required=True, help="Path to the peptides TXT file (required)")
    parser.add_argument("-m", "--mapping", type=Path, help="Optional precomputed peptide-to-protein mapping file")
    parser.add_argument("-s", "--scores", type=Path, help="Optional precomputed external scores CSV file with columns: 'Protein', 'external_score'")
    parser.add_argument("-c", "--compute_psauron", action="store_true", help="Compute PSAURON scores")
    parser.add_argument("-rg", "--reference_gtf", type=Path, help="Optional reference gtf file")
    parser.add_argument("-rf", "--reference_fasta", type=Path, help="Optional reference fasta file")
    parser.add_argument("-o", "--output", type=Path, help="Optional output directory")

    args = parser.parse_args()
    

    # Validate required arguments
    if not args.proteins:
        if not args.assembly:
            parser.error("If --proteins is not provided, --assembly must be specified.")

    try:
        if args.output:
            output_dir = args.output
            output_dir.mkdir(parents=True, exist_ok=True)
        else:
            base_dir = args.gtf.parent.parent
            output_dir = base_dir / "GAPMS_Output"
            output_dir.mkdir(parents=True, exist_ok=True)
            
        figures_dir= output_dir / "Figures"
        figures_dir.mkdir(parents=True, exist_ok=True)

        txt_dir= output_dir / "Txt"
        txt_dir.mkdir(parents=True, exist_ok=True)

        new_proteins_dir = output_dir / "Novel"
        new_proteins_dir.mkdir(parents=True, exist_ok=True)


        # Prepare logging
        log_file_path = output_dir / "log.txt"
        log_file = open(log_file_path, "w")
        sys.stdout = log_file
        sys.stderr = log_file

        # Print to terminal
        print(f"\nStarting GAP-MS. Check the log file in {log_file_path}", file=sys.__stdout__)

        # Start the pipeline
        if args.proteins:
            protein_fasta = args.proteins
        else:
            print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Generating protein FASTA from GTF and genome assembly using Gffread...")
            protein_fasta = run_gffread(args.assembly, args.gtf, output_dir, "predictions")

        if args.mapping:
            mapping = args.mapping
            print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Using provided mapping file: {mapping}")
        else:
            print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Mapping peptides to proteins with Proteomapper....")
            mapping = run_proteomapper(protein_fasta, args.peptides, output_dir)

        # Handle external scores
        external_scores_csv = None

        if args.scores:
            external_scores_csv = args.scores
            print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Using provided external scores file: {external_scores_csv}")
        elif args.compute_psauron:
            print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Computing PSAURON scores...")
            external_scores_csv = run_psauron(protein_fasta, output_dir)

        print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Extracting Features....")
        extract_features(args.gtf, protein_fasta, mapping, output_dir, external_scores_csv)

        print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Running GAPMS....")
        filter_predictions(args.gtf, protein_fasta, output_dir)

        print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Mapping peptides to genome coordinates....")
        map_peptides_to_genome(args.gtf, protein_fasta, mapping, output_dir)

        if args.reference_fasta:
            print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Getting New Supported Proteins....")
            get_new_proteins(output_dir, args.reference_gtf, args.reference_fasta)

        elif args.reference_gtf and args.assembly:
            print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Getting New Supported Proteins....")
            run_gffread(args.assembly, args.reference_gtf, output_dir, "reference")
            get_new_proteins(output_dir, args.reference_gtf)
            

        if args.reference_gtf:
            print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Getting New Supported Proteins....")
            gtf_out = output_dir / "supported_proteins.gtf"
            run_gffcompare(args.reference_gtf, gtf_out)
            parse_gffcompare_tmap(output_dir, args.gtf, args.reference_gtf, args.reference_fasta)

        total_time = (time.time() - start_time) / 60
        print("\n✅ Pipeline completed", file=sys.__stdout__)
        print(f"\nTotal time run: {total_time:.2f} minutes", file=sys.__stdout__)

    except Exception as e:
        print(f"\n❌ An error occurred: {str(e)}", file=sys.__stdout__)
        print(f"Please check the log file for more details: {log_file_path}", file=sys.__stdout__)
        
        traceback.print_exc(file=log_file)
        
        raise


    finally:
        try:
            log_file.close()
        except:
            pass

if __name__ == "__main__":
    main()
