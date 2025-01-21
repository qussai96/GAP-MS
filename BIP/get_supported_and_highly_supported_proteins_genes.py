#!/usr/bin/env python

import sys
import pandas as pd
import re
from pathlib import Path

output_dir = Path.cwd()/ 'outputs'
output_dir.mkdir(exist_ok=True)

def extract_id(suppl_value, pattern):
    match = re.search(pattern, suppl_value)
    if match:
        return match.group(1)
    return None

def gtf_to_df_with_genes(gtf_file):
    column_names = ['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Suppl']
    df = pd.read_csv(gtf_file, sep='\t', index_col=False, names=column_names, dtype={'Start': int, 'End': int})
    
    df['Gene'] = df.apply(lambda row: extract_id(row['Suppl'], r'gene_id "([^"]+)"') , axis=1)
    df['Protein'] = df.apply(lambda row: extract_id(row['Suppl'], r'transcript_id "([^"]+)"'), axis=1)
    return df

def mapping_file_to_df(mapped_file):
    mapping_df= pd.read_csv(mapped_file, sep='\t', index_col=False, names= ['Peptide', 'Protein', 'Uniqueness'])
    print(f'Total peptides loaded: {len(mapping_df)}')
    return mapping_df

def get_gene_protein_specific_peps(pred_df, mapping_df):
    # Get only list of proteins and genes from the prediction
    pred_df = pred_df[['Protein', 'Gene']]

    # Merge with mapping df to get peptide-gene df
    pep_protein_gene_df = mapping_df.merge(pred_df, on='Protein', how='left').drop_duplicates().dropna()

    # Create a set of unique proteins and genes for each peptide
    peptide_protein_map = pep_protein_gene_df.groupby('Peptide')['Protein'].nunique()
    peptide_gene_map = pep_protein_gene_df.groupby('Peptide')['Gene'].nunique()

    # Map the results back to the original dataframe
    pep_protein_gene_df['Protein_specific'] = pep_protein_gene_df['Peptide'].map(peptide_protein_map).eq(1).map({True: '+', False: '-'})
    pep_protein_gene_df['Gene_specific'] = pep_protein_gene_df['Peptide'].map(peptide_gene_map).eq(1).map({True: '+', False: '-'})

    pep_protein_gene_df.to_csv(output_dir / 'all_pep_protein_gene_df.tsv', sep='\t', index=False, header=True)
    return pep_protein_gene_df

def get_supported_genes_proteins(pep_protein_gene_df):
    # Filter and save supported proteins and genes
    supported_proteins_df = pep_protein_gene_df[pep_protein_gene_df['Protein_specific'] == '+']
    supported_proteins_df.to_csv(output_dir / 'supported_proteins.tsv', sep='\t', index=False, header=True)
    print(f'Number of supported proteins: {len(set(supported_proteins_df["Protein"]))}')

    supported_genes_df = pep_protein_gene_df[pep_protein_gene_df['Gene_specific'] == '+']
    supported_genes_df.to_csv(output_dir / 'supported_genes.tsv', sep='\t', index=False, header=True)
    print(f'Number of supported genes: {len(set(supported_genes_df["Gene"]))}')


def get_highly_supported_genes_proteins(pep_protein_gene_df):
    # Count occurrences of protein_specific and gene_specific per protein and gene
    protein_specific_count = pep_protein_gene_df.groupby('Protein')['Protein_specific'].apply(lambda x: (x == '+').sum())
    gene_specific_count = pep_protein_gene_df.groupby('Gene')['Gene_specific'].apply(lambda x: (x == '+').sum())

    # Filter proteins and genes with at least two protein_specific or gene_specific peptides
    highly_supported_proteins_df = pep_protein_gene_df[pep_protein_gene_df['Protein'].isin(protein_specific_count[protein_specific_count >= 2].index)]
    highly_supported_genes_df = pep_protein_gene_df[pep_protein_gene_df['Gene'].isin(gene_specific_count[gene_specific_count >= 2].index)]

    # Output the DataFrames to CSV
    highly_supported_proteins_df.to_csv(output_dir / 'highly_supported_proteins.tsv', sep='\t', index=False, header=True)
    highly_supported_genes_df.to_csv(output_dir / 'highly_supported_genes.tsv', sep='\t', index=False, header=True)

    print(f'Number of highly supported proteins: {len(set(highly_supported_proteins_df["Protein"]))}')
    print(f'Number of highly supported genes: {len(set(highly_supported_genes_df["Gene"]))}')






def main(argv):
    pred_gtf = argv[0]
    mapping_file = argv[1]
    
    pred_df = gtf_to_df_with_genes(pred_gtf)
    mapping_df = mapping_file_to_df(mapping_file)
    
    pep_protein_gene_df = get_gene_protein_specific_peps(pred_df, mapping_df)
    get_supported_genes_proteins(pep_protein_gene_df)
    get_highly_supported_genes_proteins(pep_protein_gene_df)

if __name__ == "__main__":
    main(sys.argv[1:])
