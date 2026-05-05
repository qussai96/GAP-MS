import numpy as np
import pandas as pd

def mapping_file_to_df(mapping_file):
    mapping_df = pd.read_csv(mapping_file, sep='\t', header=0, names=['Peptide', 'Protein', 'location'])
    mapping_df = mapping_df[['Peptide', 'Protein']]
    mapping_df.drop_duplicates(inplace=True)
    mapping_df = mapping_df[~mapping_df['Protein'].str.lower().eq('unmapped')]
    mapping_df['peptide_length'] = mapping_df['Peptide'].apply(len)
    print(f'Total peptides loaded: {mapping_df["Peptide"].nunique()}')
    return mapping_df


def get_gene_protein_specific_peps(gtf_df, mapping_df):
    pred_df = gtf_df[['Protein', 'Gene']]
    pep_protein_gene_df = mapping_df.merge(pred_df, on='Protein', how='left').drop_duplicates().dropna()
    
    peptide_protein_map = pep_protein_gene_df.groupby('Peptide')['Protein'].nunique()
    peptide_gene_map = pep_protein_gene_df.groupby('Peptide')['Gene'].nunique()
    pep_protein_gene_df['Protein_specific'] = pep_protein_gene_df['Peptide'].map(peptide_protein_map).eq(1).map({True: '+', False: '-'})
    pep_protein_gene_df['Gene_specific'] = pep_protein_gene_df['Peptide'].map(peptide_gene_map).eq(1).map({True: '+', False: '-'})
    return pep_protein_gene_df[['Peptide', 'Protein', 'Gene', 'peptide_length', 'Protein_specific', 'Gene_specific']]



def calculate_sequence_coverage(peptide_df):
    """
    Compute peptide-derived sequence coverage per protein.

    Accepts either a peptide-level DataFrame containing multiple proteins and returns
    a `{protein_id: coverage}` dictionary, or a single-protein grouped DataFrame and
    returns the coverage as a float.
    """
    required_cols = {'Protein', 'pep_start', 'pep_end', 'prot_len'}
    missing_cols = required_cols - set(peptide_df.columns)
    if missing_cols:
        raise ValueError(f"Missing required columns for coverage calculation: {sorted(missing_cols)}")

    if peptide_df.empty:
        return {}

    valid_peptides = peptide_df.loc[
        peptide_df['pep_start'].fillna(0).gt(0)
        & peptide_df['pep_end'].fillna(0).ge(peptide_df['pep_start'].fillna(0))
        & peptide_df['prot_len'].fillna(0).gt(0),
        ['Protein', 'pep_start', 'pep_end', 'prot_len']
    ].copy()

    all_proteins = peptide_df['Protein'].dropna().unique().tolist()
    if valid_peptides.empty:
        coverage_map = {protein: 0.0 for protein in all_proteins}
        return coverage_map[all_proteins[0]] if len(all_proteins) == 1 else coverage_map

    pep_sorted = valid_peptides.sort_values(['Protein', 'pep_start', 'pep_end']).reset_index(drop=True)

    prot_arr = pep_sorted['Protein'].to_numpy()
    start_arr = pep_sorted['pep_start'].to_numpy()
    end_arr = pep_sorted['pep_end'].to_numpy()
    len_arr = pep_sorted['prot_len'].to_numpy()

    cur_prot = None
    cur_start = cur_end = 0
    total = 0
    prot_len_val = 0
    coverage_map = {}

    for i in range(len(prot_arr)):
        p = prot_arr[i]
        s = start_arr[i]
        e = end_arr[i]

        if p != cur_prot:
            if cur_prot is not None:
                total += cur_end - cur_start + 1
                coverage_map[cur_prot] = round(total / prot_len_val, 4)
            cur_prot = p
            cur_start = s
            cur_end = e
            total = 0
            prot_len_val = len_arr[i]
        elif s <= cur_end:
            cur_end = max(cur_end, e)
        else:
            total += cur_end - cur_start + 1
            cur_start = s
            cur_end = e

    if cur_prot is not None:
        total += cur_end - cur_start + 1
        coverage_map[cur_prot] = round(total / prot_len_val, 4)

    for protein in all_proteins:
        coverage_map.setdefault(protein, 0.0)

    return coverage_map[all_proteins[0]] if len(all_proteins) == 1 else coverage_map


