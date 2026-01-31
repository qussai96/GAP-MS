import sys
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from .gtf_utils import gtf_to_df_with_genes


def parse_gffcompare_tmap(output_dir, prediction_gtf, reference_gtf, reference_fasta=None):
    # 1. Setup Paths
    output_dir = Path(output_dir)
    new_proteins_dir = output_dir / "Novel"
    new_proteins_dir.mkdir(parents=True, exist_ok=True)

    # Set file paths
    # Note: Adjust tmap filename if your gffcompare output prefix is different
    # This looks for any .tmap file in the dir, or defaults to specific name
    tmap_files = list(output_dir.glob("*.tmap"))
    if tmap_files:
        tmap_file = tmap_files[0]
    else:
        tmap_file = output_dir / "gffcompare.supported_proteins.gtf.tmap"
    
    supported_proteins_fasta = output_dir / "supported_proteins.fa"
    all_proteins_scores = output_dir / "all_proteins_scores.tsv"
    
    reference_fasta = Path(reference_fasta) if reference_fasta else output_dir / "translated_reference_proteins.fa"

    # 2. Parse gffcompare TMAP file
    print(f"Parsing TMAP file: {tmap_file}")
    # TMAP is tab-separated. Key columns: 
    # 'class_code' (col 3) and 'qry_id' (col 5, which matches your Protein ID)
    tmap_df = pd.read_csv(tmap_file, sep='\t')
    
    # Filter for Identical proteins (Class code '=')
    identical_proteins = set(tmap_df[tmap_df['class_code'] == '=']['qry_id'])
    
    # Filter for New/Novel proteins based on your requested codes
    # u: unknown/intergenic, x: antisense, i: intronic, etc.
    novel_codes = ['s', 'x', 'i', 'y', 'p', 'u']
    new_predicted_df = tmap_df[tmap_df['class_code'].isin(novel_codes)]
    new_predicted_proteins = set(new_predicted_df['qry_id'])

    # 3. Print Statistics
    # We count FASTA records for accurate file lengths
    ref_count = len(list(SeqIO.parse(reference_fasta, "fasta")))
    supp_count = len(list(SeqIO.parse(supported_proteins_fasta, "fasta")))

    print(f'Number of reference transcripts = {ref_count}')
    print(f'Number of supported predicted transcripts = {supp_count}')
    print(f'Number of identical predicted transcripts = {len(identical_proteins)}')
    print(f'Number of newly predicted transcripts = {len(new_predicted_proteins)}')

    # 4. Generate GTF for New Proteins
    print("Generating GTF for new proteins...")
    prediction_gtf_df = gtf_to_df_with_genes(prediction_gtf)
    # Filter for CDS only
    prediction_gtf_df = prediction_gtf_df[prediction_gtf_df['Type'].isin(['CDS'])]
    
    # Keep only rows belonging to new proteins
    new_predicted_proteins_gtf = prediction_gtf_df[prediction_gtf_df['Protein'].isin(new_predicted_proteins)]
    
    # Remove 'Protein' and 'Gene' columns before saving (as per your snippet)
    cols_to_drop = [c for c in ['Protein', 'Gene'] if c in new_predicted_proteins_gtf.columns]
    new_predicted_proteins_gtf = new_predicted_proteins_gtf.drop(cols_to_drop, axis=1)
    
    new_predicted_proteins_gtf.to_csv(new_proteins_dir / "new_predicted_proteins.gtf", sep='\t', header=False, index=False)

    # 5. Generate FASTA for New Proteins
    print("Generating FASTA for new proteins...")
    prediction_seq_dict = {record.id: str(record.seq) for record in SeqIO.parse(supported_proteins_fasta, "fasta")}
    
    # Extract sequences for new proteins only
    new_predicted_proteins_dict = {
        key: prediction_seq_dict[key] 
        for key in new_predicted_proteins 
        if key in prediction_seq_dict
    }

    def write_seq_dict_to_file(data, filename):
        with open(filename, "w") as file:
            for key, value in data.items():
                file.write(f">{key}\n{value}\n")
                
    write_seq_dict_to_file(new_predicted_proteins_dict, new_proteins_dir / "new_predicted_proteins.fa")

    # 6. Extract Scores for New Proteins
    print("Extracting scores for new proteins...")
    if all_proteins_scores.exists():
        all_protein_features_df = pd.read_csv(all_proteins_scores, sep='\t', header=0)
        
        # Ensure the score file has a 'Protein' column to match against
        if 'Protein' in all_protein_features_df.columns:
            new_predicted_proteins_df = all_protein_features_df[
                all_protein_features_df['Protein'].isin(new_predicted_proteins)
            ]
            new_predicted_proteins_df.to_csv(new_proteins_dir / "new_predicted_proteins_scores.tsv", sep='\t', index=False)
        else:
            print("Warning: 'Protein' column missing in scores file. Skipping score extraction.")
    else:
        print(f"Warning: Scores file not found at {all_proteins_scores}")

    # 7. Create Combined FASTA (Reference + New)
    print("Creating combined FASTA (Reference + New)...")
    sequences = []
    
    # Load Reference
    sequences.extend(SeqIO.parse(reference_fasta, "fasta"))
    
    # Load New Proteins (from the file we just created)
    sequences.extend(SeqIO.parse(new_proteins_dir / "new_predicted_proteins.fa", "fasta"))

    # Write combined
    combined_fasta_path = new_proteins_dir / "reference+new_proteins.fa"
    with open(combined_fasta_path, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

    print(f"Combined FASTA file saved at: {combined_fasta_path}")