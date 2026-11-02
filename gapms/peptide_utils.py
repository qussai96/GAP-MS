import numpy as np
import pandas as pd


def resolve_peptide_starts(protein_sequence, peptide, location=None):
    """Return validated 1-based starts for a peptide-protein mapping.

    The built-in protein mapper reports one-based positions.  A zero-based compatibility check
    is retained for user-supplied mapping files, but one-based wins when both
    interpretations happen to match (for example in a homopolymer).  If no
    usable position is supplied, return every sequence occurrence rather than
    silently relocating a repeated peptide to its first occurrence.
    """
    sequence = str(protein_sequence).replace('I', 'L').replace('*', '')
    normalized_peptide = str(peptide).replace('I', 'L')
    if not normalized_peptide:
        return []

    try:
        numeric_location = int(location)
    except (TypeError, ValueError, OverflowError):
        numeric_location = None

    if numeric_location is not None:
        if (
            numeric_location >= 1
            and sequence[
                numeric_location - 1:numeric_location - 1 + len(normalized_peptide)
            ] == normalized_peptide
        ):
            return [numeric_location]
        if (
            numeric_location >= 0
            and sequence[numeric_location:numeric_location + len(normalized_peptide)]
            == normalized_peptide
        ):
            return [numeric_location + 1]

    starts = []
    cursor = 0
    while True:
        index = sequence.find(normalized_peptide, cursor)
        if index < 0:
            break
        starts.append(index + 1)
        cursor = index + 1
    return starts

def mapping_file_to_df(mapping_file):
    mapping_df = pd.read_csv(mapping_file, sep='\t', header=0, names=['Peptide', 'Protein', 'location'])
    mapping_df = mapping_df[['Peptide', 'Protein', 'location']]
    mapping_df.drop_duplicates(inplace=True)
    mapping_df = mapping_df[~mapping_df['Protein'].str.lower().eq('unmapped')]
    
    # Split multiple proteins on same line (separated by spaces) into separate rows
    mapping_df['Protein'] = mapping_df['Protein'].str.split()
    mapping_df = mapping_df.explode('Protein')
    mapping_df['location'] = pd.to_numeric(mapping_df['location'], errors='coerce')
    
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
    columns = [
        'Peptide', 'Protein', 'Gene', 'peptide_length',
        'Protein_specific', 'Gene_specific',
    ]
    if 'location' in pep_protein_gene_df.columns:
        columns.append('location')
    return pep_protein_gene_df[columns]



def calculate_sequence_coverage(peptide_df):
    """
    Compute peptide-derived sequence coverage per protein.

    Returns a `{protein_id: coverage}` dictionary, including for a one-protein
    input so callers do not need a dataset-size-dependent special case.
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
        return coverage_map

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

    return coverage_map
