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


def check_peptide_loc(row, protein_cds_dict):
    protein = row['Protein']
    pep_start, pep_end = row['pep_start'], row['pep_end']
    protein_cds = protein_cds_dict.get(protein, [])
    
    for cds_start, cds_end in protein_cds:
        if cds_start <= pep_start and cds_end >= pep_end:
            return '-'
    return '+'

def calculate_protein_coverage(group):
    prot_len = group["prot_len"].iloc[0]
    intervals = group[["pep_start", "pep_end"]].to_numpy()
    intervals = intervals[np.argsort(intervals[:, 0])]
    
    merged_intervals = []
    start, end = intervals[0]

    for s, e in intervals[1:]:
        if s <= end:
            end = max(end, e)
        else:
            merged_intervals.append((start, end))
            start, end = s, e
    merged_intervals.append((start, end))

    covered_length = sum(e - s + 1 for s, e in merged_intervals)
    return round(covered_length / prot_len, 4)

def count_expected_peptides_with_missed_cleavages(sequence, missed_cleavages=1):
    cleavage_sites = [
        i for i in range(len(sequence) - 1)
        if sequence[i] in ['K', 'R'] and sequence[i + 1] != 'P'
    ]
    num_cleavages = len(cleavage_sites)
    total_peptides = sum((num_cleavages + missed) for missed in range(missed_cleavages + 1))
    return total_peptides

