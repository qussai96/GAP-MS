#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unites two gene predictions
- It should be specified which one is the main prediction.
  All entries from the main prediction will have a corresponding label, but only
  new proteins from the additional prediction that are not present in the main
  prediction will have a label of the additional prediction.
"""

import pandas as pd
import re
import os
import sys
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
    print(f'Predicted proteins loaded from {gtf_file}: {len(set(df["Protein"]))}')
    return df

def mapping_file_to_df(mapped_file):
    mapping_df= pd.read_csv(mapped_file, sep='\t', index_col=False, names= ['Peptide', 'Protein', 'Protein_specific'])
    print(f'Total peptides loaded: {len(mapping_df)}')
    return mapping_df


def len_overlap(cds1, cds2):
    """
    Calculate the nucleotide overlap between two sets of CDS intervals.

    Args:
        cds1 (list of lists): CDS intervals for the first gene, e.g., [[start1, end1], [start2, end2], ...].
        cds2 (list of lists): CDS intervals for the second gene, in the same format.

    Returns:
        int: Total overlap in nucleotides between the two sets of CDS intervals.
    """
    # Convert to DataFrame for vectorized operations
    cds1_df = pd.DataFrame(cds1, columns=['Start', 'End'])
    cds2_df = pd.DataFrame(cds2, columns=['Start', 'End'])

    # Create all combinations of intervals and calculate overlaps
    overlaps = pd.merge(cds1_df.assign(key=1), cds2_df.assign(key=1), on='key')
    overlaps['OverlapStart'] = overlaps[['Start_x', 'Start_y']].max(axis=1)
    overlaps['OverlapEnd'] = overlaps[['End_x', 'End_y']].min(axis=1)
    overlaps = overlaps[overlaps['OverlapStart'] <= overlaps['OverlapEnd']]

    # Calculate total overlap
    total_overlap = (overlaps['OverlapEnd'] - overlaps['OverlapStart'] + 1).sum()
    return total_overlap


def rename_genes_proteins(df, prefix):
    """Add a prefix to gene and protein names to avoid conflicts."""
    df['Protein'] = prefix + df['Protein']
    df['Gene'] = prefix + df['Gene']
    return df

def find_overlapping_genes(df1, df2, min_overlap):
    overlapping_genes1, overlapping_genes2, new_genes2 = [], [], []

    # Group df2 by gene for efficient processing
    df2_grouped = df2.groupby('Gene')

    for gene, df2_gene in df2_grouped:
        seqid, strand = df2_gene.iloc[0][['Seqid', 'Strand']]
        start, end = df2_gene['Start'].min(), df2_gene['End'].max()

        # Filter df1 for potential overlaps
        df1_candidates = df1[
            (df1['Seqid'] == seqid) &
            (df1['Strand'] == strand) &
            (df1['Start'] < end) &
            (df1['End'] > start)
        ]

        if df1_candidates.empty:
            new_genes2.append(gene)
            continue

        # Calculate overlaps for all candidates
        df1_grouped = df1_candidates.groupby('Gene')
        overlaps = {
            candidate_gene: len_overlap(
                df1_grouped.get_group(candidate_gene)[['Start', 'End']].values.tolist(),
                df2_gene[['Start', 'End']].values.tolist()
            )
            for candidate_gene in df1_grouped.groups.keys()
        }

        # Find the best overlap
        overlaps = {gene: overlap for gene, overlap in overlaps.items() if overlap >= min_overlap}
        if overlaps:
            best_gene = max(overlaps, key=overlaps.get)
            overlapping_genes1.append(best_gene)
            overlapping_genes2.append(gene)
        else:
            new_genes2.append(gene)

    return overlapping_genes1, overlapping_genes2, new_genes2

def update_gene_names(df2, overlapping_genes1, overlapping_genes2):
    """Update gene names in df2 based on overlapping genes in df1."""
    name_mapping = dict(zip(overlapping_genes2, overlapping_genes1))
    df2['Gene'] = df2['Gene'].apply(lambda x: name_mapping.get(x, x))
    return df2

def filter_unique_proteins(df1, df2):
    print("filter")
    """Remove proteins in df2 that are already present in df1."""
    df1_proteins = set(
        df1.groupby('Protein').apply(lambda x: f"{x.iloc[0]['Seqid']}:{','.join(map(str, x[['Start', 'End']].values.flatten()))}")
    )

    df2['Protein_info'] = df2.groupby('Protein').apply(lambda x: f"{x.iloc[0]['Seqid']}:{','.join(map(str, x[['Start', 'End']].values.flatten()))}")
    df2 = df2[~df2['Protein_info'].isin(df1_proteins)]
    df2 = df2.drop(columns=['Protein_info'])
    return df2



def unite_predictions(df1, df2, min_overlap, pred_name1, pred_name2):
    """Combine two prediction dataframes, resolving overlaps and unique entries."""
    # Step 1: Rename genes and proteins in df2 to avoid conflicts
    print("Rename genes and proteins....")
    df2 = rename_genes_proteins(df2, 'x')

    # Step 2: Find overlapping and new genes
    print("Finding overlapping genes....")
    overlapping_genes1, overlapping_genes2, new_genes2 = find_overlapping_genes(df1, df2, min_overlap)

    # Step 3: Update gene names in df2 for overlapping genes
    print("Update gene names....")
    df2 = update_gene_names(df2, overlapping_genes1, overlapping_genes2)

    # Step 4: Filter out proteins in df2 already present in df1
    print("Filter out proteins....")
    df2 = filter_unique_proteins(df1, df2)

    # Step 5: Combine df1 and df2
    print("Combine predictions....")
    df1['Prediction'] = pred_name1
    df2['Prediction'] = pred_name2
    print("step5")
    combined_df = pd.concat([df1, df2], ignore_index=True)

    # Step 6: Update the 'Suppl' column and sort
    combined_df['Suppl'] = combined_df.apply(
        lambda row: f"transcript_id \"{row['Protein']}\"; gene_id \"{row['Gene']}\"; prediction \"{row['Prediction']}\";",
        axis=1
    )
    combined_df = combined_df[['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Suppl']]
    combined_df = combined_df.sort_values(['Seqid', 'Start', 'End']).reset_index(drop=True)
    combined_df.to_csv(output_dir / f'final_predictions.gtf', sep='\t', index=False, header=False) 



def main(argv):
    file_main_pred = argv[0]
    file_second_pred = argv[1]
    label_main_pred = argv[2]
    label_second_pred = argv[3]
    min_overlap = int(argv[4])

    main_pred_df = gtf_to_df_with_genes(file_main_pred)
    second_pred_df = gtf_to_df_with_genes(file_second_pred)
    unite_predictions(main_pred_df, second_pred_df, min_overlap, label_main_pred, label_second_pred)


if __name__ == "__main__":
    main(sys.argv[1:])
