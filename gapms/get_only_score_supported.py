import re
import pandas as pd
import csv
import sys
from pathlib import Path



def extract_id(suppl_value, pattern):
    if isinstance(suppl_value, str):
        match = re.search(pattern, suppl_value)
        if match:
            return match.group(1)
    # print(f'Check row with suppl {suppl_value}')
    return None

def gtf_to_df_with_genes(gtf_file):
    column_names = [
        'Seqid', 'Source', 'Type', 'Start', 'End',
        'Prediction_score', 'Strand', 'Frame', 'Suppl'
    ]

    # Read GTF/GFF file
    df = pd.read_csv(
        gtf_file,
        sep='\t',
        header=None,         # prevent first line being treated as header
        names=column_names,
        comment='#'
    )

    df['Type'] = df['Type'].astype(str)
    df = df[df['Type'].str.lower().isin(['gene', 'mrna', 'transcript', 'cds'])]
    def gene_extraction(row):
        if "gene_id" in row['Suppl']:
            return extract_id(row['Suppl'], r'gene_id "([^"]+)"')
        elif "Parent" in row['Suppl']:
            return extract_id(row['Suppl'], r'Parent=([^;]+)')
        elif row['Type'] == 'gene' and "ID" in row['Suppl']:
            return extract_id(row['Suppl'], r'ID=([^;]+)')
        elif row['Type'].lower() in ['mrna', 'transcript', 'gene']:
            return row['Suppl']
        else:
            print(f'No gene id found. Check row with suppl {row}')
            return None

    def protein_extraction(row):
        if row['Type'].lower() in ['mrna', 'transcript', 'gene']:
            if "ID=" in row['Suppl']:
                return extract_id(row['Suppl'], r'ID=([^;]+)')
            elif "transcript_id" in row['Suppl']:
                return extract_id(row['Suppl'], r'transcript_id "([^"]+)"')
            elif row['Type'] == 'transcript':
                return row['Suppl']
            else:
                return row['Suppl']
        else:
            if "protein_id" in row['Suppl']:
                return extract_id(row['Suppl'], r'protein_id "([^"]+)"')
            elif "transcript_id" in row['Suppl']:
                return extract_id(row['Suppl'], r'transcript_id "([^"]+)"')
            elif "Parent" in row['Suppl']:
                return extract_id(row['Suppl'], r'Parent=([^;]+)')
            else:
                print(f'No protein id found. Check row with suppl {row}')
                return None
    
    df = df.dropna(subset=['Seqid'])
    df['Gene'] = df.apply(gene_extraction, axis=1)
    df['Protein'] = df.apply(protein_extraction, axis=1)

    # Convert prediction scores safely
    df['Prediction_score'] = pd.to_numeric(
        df['Prediction_score'].replace('.', None, regex=False),
        errors='coerce'
    ).fillna(0)
    return df



def save_gtf_subset(gtf_df, subset, output_path):
    output_path = Path(output_path)
    subset_gtf = gtf_df[gtf_df['Protein'].isin(subset)]
    subset_gtf.drop(['Protein', 'Gene'], axis=1).to_csv(
        output_path,
        sep="\t",
        index=False,
        header=False,
        quoting=csv.QUOTE_NONE,
        escapechar=' '
    )

def filter_for_score_supported(gtf, scores_tsv, output_path, score_threshold=0.5):
    # Convert GTF to DataFrame
    gtf_df = gtf_to_df_with_genes(gtf)
    
    # Load scores table
    scores_df = pd.read_csv(
        scores_tsv, sep=",", names=["Protein", "is_protein", "score"],
        on_bad_lines="skip",
        skip_blank_lines=True).drop(columns="is_protein")
    scores_df["score"] = pd.to_numeric(scores_df["score"], errors="coerce")
    # Logging counts
    print(f"Number of all proteins = {(scores_df['Protein'].nunique())}")
    print(f"Number of supported proteins = {(scores_df[scores_df['score'] >= score_threshold]['Protein'].nunique())}")
    
    # Get supported proteins
    supported_proteins = set(scores_df[scores_df['score'] >= score_threshold]["Protein"])
    
    # Save subs
    save_gtf_subset(gtf_df, supported_proteins, output_path)

def main(argv):
    gtf = argv[0]
    scores_tsv = argv[1]
    output_path = argv[2]
    filter_for_score_supported(gtf, scores_tsv, output_path)

if __name__ == "__main__":
    main(sys.argv[1:])
