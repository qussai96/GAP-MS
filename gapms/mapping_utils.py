import pandas as pd

def mapping_file_to_df(mapping_file):
    """
    Loads the peptide-to-protein mapping file and filters out unmapped entries.

    Parameters:
    - mapping_file (str or Path): Path to the mapping file with columns [Peptide, Protein, Location]

    Returns:
    - DataFrame with columns [Peptide, Protein, Location, pep_len]
    """
    mapping_df = pd.read_csv(mapping_file, sep='\t', header=None, names=['Peptide', 'Protein', 'location'])
    mapping_df = mapping_df[~mapping_df['Protein'].str.lower().eq('unmapped')]
    mapping_df['pep_len'] = mapping_df['Peptide'].apply(len)
    print(f'Total peptides loaded: {mapping_df["Peptide"].nunique()}')
    return mapping_df
