#!/usr/bin/env python

"""
Finds a BRAKER transcript Score cutoff based on the Scores of supported Proteins, 
applies the Score filter, and makes a .gtf file with high-scoring predictions 
and all supported Proteins. 
"""

import sys
import pandas as pd
import numpy as np
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



def find_cutoff(Scores_df, sp_Proteins):
    """
    Calculate the optimal BRAKER Score cutoff based on the ROC curve.

    Parameters:
    - Scores_df (DataFrame): A DataFrame containing Protein IDs and BRAKER Scores.
    - sp_Proteins (list): A list of supported Protein IDs.

    Returns:
    - float: The optimal BRAKER Score cutoff.
    """
    # Ensure inputs are valid
    if Scores_df.empty or not sp_Proteins:
        raise ValueError("Input DataFrame is empty or supported Proteins list is empty.")

    # Extract Scores of all Proteins
    Scores = Scores_df['Score'].values.tolist()

    # Extract Scores of supported Proteins
    sp_Scores_df = Scores_df[Scores_df['Protein'].isin(sp_Proteins)]
    sp_Scores = sp_Scores_df['Score'].values.tolist()

    if not sp_Scores:
        raise ValueError("No supported Protein Scores found in the DataFrame.")

    # Generate possible cutoff Scores: from 1.0 to 0.0 with a step of 0.01
    cutoffs = [x / 100 for x in range(100, -1, -1)]

    # Calculate sensitivity (TPR) and specificity (TNR) for each cutoff
    tpr = [calculate_sensitivity(sp_Scores, cutoff) for cutoff in cutoffs]
    tnr = [calculate_specificity(Scores, sp_Scores, cutoff) for cutoff in cutoffs]
    fpr = [1 - x for x in tnr]  # False Positive Rate

    # Find the cutoff where the slope of the ROC curve drops below 1
    for i in range(2, len(tpr)):
        slope = (tpr[i] - tpr[i - 1]) / (fpr[i] - fpr[i - 1]) if (fpr[i] - fpr[i - 1]) != 0 else float('inf')
        if slope < 1:
            cutoff = cutoffs[i]
            tpr_cut = tpr[i]
            fpr_cut = fpr[i]
            break
    else:
        raise ValueError("No valid cutoff found where the slope of the ROC curve drops below 1.")

    print(f"BRAKER Score cutoff: {cutoff}")
    print(f"TPR at cutoff: {tpr_cut}")
    print(f"FPR at cutoff: {fpr_cut}")

    return cutoff

def calculate_sensitivity(Scores, cutoff):
    """
    Calculate sensitivity (True Positive Rate) at a given cutoff.

    Parameters:
    - Scores (list): List of Scores of supported Proteins.
    - cutoff (float): The cutoff Score.

    Returns:
    - float: Sensitivity as the fraction of supported Proteins with Scores >= cutoff.
    """
    if not Scores:
        return 0.0

    passed = [Score for Score in Scores if Score >= cutoff]
    return len(passed) / len(Scores)

def calculate_specificity(all_Scores, sp_Scores, cutoff):
    """
    Calculate specificity (True Negative Rate) at a given cutoff.

    Parameters:
    - all_Scores (list): List of Scores of all Proteins.
    - sp_Scores (list): List of Scores of supported Proteins.
    - cutoff (float): The cutoff Score.

    Returns:
    - float: Specificity as the fraction of supported Proteins with Scores >= cutoff among all Proteins with Scores >= cutoff.
    """
    passed_all = [Score for Score in all_Scores if Score >= cutoff]
    passed_sp = [Score for Score in sp_Scores if Score >= cutoff]

    if not passed_all:
        return 0.0

    return len(passed_sp) / len(passed_all)


def find_apply_Score_filter(gtf_df, supp_prot_df):
    """
    Apply a Score filter to identify high-scoring and supported Proteins, and generate a list of final Proteins.

    Parameters:
    - Scores_file (str): Path to the file containing Protein Scores.
    - supp_prot_df (DataFrame): DataFrame containing supported Proteins.
    - gtf_df (DataFrame): DataFrame containing GTF data (not used in this function).

    Returns:
    - list: Final list of Proteins passing the Score filter.
    """

    # Extract supported Proteins
    if 'Protein' not in supp_prot_df.columns:
        raise ValueError("The supported Proteins DataFrame must contain a 'Protein' column.")
    sp_Proteins = set(supp_prot_df['Protein'].dropna().values)

    if not sp_Proteins:
        raise ValueError("No supported Proteins found in the provided DataFrame.")

    
    # Replace missing or invalid Scores with a default value and convert to float
    Scores_df = gtf_df[gtf_df['Type'] == 'transcript']
    Scores_df = Scores_df[['Suppl', 'Score']]
    Scores_df.rename(columns={"Suppl": "Protein"}, inplace=True)

    Scores_df['Score'] = Scores_df['Score'].replace('.', '0.001').astype(float)
    
    
    # Find BRAKER Score cutoff
    print('Calculating the BRAKER transcript Score cutoff...')
    cutoff = find_cutoff(Scores_df, sp_Proteins)

    # Identify high-scoring Proteins and low-scoring supported Proteins
    df_high = Scores_df[Scores_df['Score'] >= cutoff]
    prot_high = set(df_high['Protein'].values)

    df_low = Scores_df[Scores_df['Score'] < cutoff]
    prot_low = set(df_low[df_low['Protein'].isin(sp_Proteins)]['Protein'].values)

    # Combine high-scoring Proteins and low-scoring supported Proteins
    final_Proteins = prot_high.union(prot_low)

    # Output statistics
    print(f'High-scoring Proteins: {len(prot_high)}')
    print(f'Low-scoring supported Proteins: {len(prot_low)}')
    print(f'Total final Proteins: {len(final_Proteins)}')

    return list(final_Proteins)



def main(argv):
    gtf_file = argv[0]
    supp_prot_gtf = argv[1]
    name = argv[2]
    
    gtf_df = gtf_to_df_with_genes(gtf_file)
    sp_prot_df = gtf_to_df_with_genes(supp_prot_gtf)
    final_Proteins = find_apply_Score_filter(gtf_df, sp_prot_df)


    final_df = gtf_df[gtf_df['Protein'].isin(final_Proteins)]
    final_df = final_df[['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Suppl']]
    final_df.to_csv(output_dir / f'high_scoring_{name}.gtf', sep='\t', index=False, header=False)         


if __name__ == "__main__":
    main(sys.argv[1:])
