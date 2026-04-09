import re
import pandas as pd
import csv


def extract_attribute(suppl_value, attr_name):
    """
    Extract attribute value from GFF/GTF attributes field.
    Handles both GTF format (attr_name "value") and GFF format (attr_name=value).
    """
    if not isinstance(suppl_value, str):
        return None
    
    # Try GTF format first: attr_name "value"
    match = re.search(rf'{attr_name}\s*"([^"]+)"', suppl_value)
    if match:
        return match.group(1)
    
    # Try GFF format: attr_name=value
    match = re.search(rf'{attr_name}=([^;,\s]+)', suppl_value)
    if match:
        return match.group(1)
    
    return None

def _str_extract(series, pattern):
    """Return first capture group from a Series using vectorized regex."""
    return series.str.extract(pattern, expand=False)


def gtf_to_df_with_genes(gtf_file):
    column_names = [
        'Seqid', 'Source', 'Type', 'Start', 'End',
        'Prediction_score', 'Strand', 'Frame', 'Suppl'
    ]

    # Read GTF/GFF file
    df = pd.read_csv(
        gtf_file,
        sep='\t',
        header=None,
        names=column_names,
        comment='#',
        dtype={'Type': str, 'Suppl': str},
    )

    df['Type'] = df['Type'].str.lower()
    df = df[df['Type'].isin(['gene', 'mrna', 'transcript', 'cds'])]
    df = df.dropna(subset=['Seqid'])

    s = df['Suppl']
    t = df['Type']

    # Vectorized gene extraction:
    # gene_id > GeneID > gene > Parent (non-gene) > ID (gene) > raw attributes fallback
    gene = (
        _str_extract(s, r'gene_id\s*[= ]?"?([^";,\s]+)"?')
        .combine_first(_str_extract(s, r'GeneID=([^;,\s]+)'))
        .combine_first(_str_extract(s, r'(?:^|;)\s*gene\s*[= ]?"?([^";,\s]+)"?'))
    )

    non_gene_mask = t != 'gene'
    gene = gene.combine_first(
        _str_extract(s.where(non_gene_mask), r'Parent=([^;,\s]+)')
    )

    gene_row_mask = t == 'gene'
    gene = gene.combine_first(
        _str_extract(s.where(gene_row_mask), r'ID=([^;,\s]+)')
    )

    meta_mask = t.isin(['mrna', 'transcript', 'gene'])
    gene = gene.fillna(s.where(meta_mask))
    df['Gene'] = gene

    # Vectorized protein extraction:
    # transcript/mRNA/gene rows: ID > transcript_id > raw attributes
    prot_meta = (
        _str_extract(s.where(meta_mask), r'ID=([^;,\s]+)')
        .combine_first(_str_extract(s.where(meta_mask), r'transcript_id\s*[= ]?"?([^";,\s]+)"?'))
        .fillna(s.where(meta_mask))
    )

    # CDS rows: protein_id > transcript_id > Parent
    cds_mask = t == 'cds'
    prot_cds = (
        _str_extract(s.where(cds_mask), r'protein_id\s*[= ]?"?([^";,\s]+)"?')
        .combine_first(_str_extract(s.where(cds_mask), r'transcript_id\s*[= ]?"?([^";,\s]+)"?'))
        .combine_first(_str_extract(s.where(cds_mask), r'Parent=([^;,\s]+)'))
    )
    df['Protein'] = prot_meta.combine_first(prot_cds)

    # Normalize Type values to what downstream code expects
    type_map = {'gene': 'gene', 'mrna': 'mRNA', 'transcript': 'transcript', 'cds': 'CDS'}
    df['Type'] = df['Type'].map(type_map).fillna(df['Type'])

    # Convert prediction scores safely
    df['Prediction_score'] = pd.to_numeric(
        df['Prediction_score'].replace('.', None, regex=False),
        errors='coerce'
    ).fillna(0)

    return df



def save_gtf_subset(gtf_df, subset, output_dir, output_file_name):
    subset_gtf = gtf_df[gtf_df['Protein'].isin(subset)]
    subset_gtf.drop(['Protein', 'Gene'], axis=1).to_csv(
        output_dir / output_file_name,
        sep="\t",
        index=False,
        header=False,
        quoting=csv.QUOTE_NONE,
        escapechar=' '
    )

def read_scores_csv(scores_csv):
    # Try reading the file while skipping metadata lines (no commas)
    with open(scores_csv, 'r') as f:
        lines = f.readlines()
    # Find first line that looks like a header (has a comma)
    start_idx = next((i for i, line in enumerate(lines) if ',' in line), 0)
    
    # Read CSV from that line
    df = pd.read_csv(scores_csv, skiprows=start_idx, on_bad_lines='skip', skip_blank_lines=True)

    # Standardize column names
    df.columns = df.columns.str.strip().str.lower()

    # Handle possible column name variants
    col_map = {
        'description': 'Protein',
        'Protein': 'Protein',
        'external_score': 'external_score',
        'in-frame_score': 'external_score'
    }

    df = df.rename(columns={old: col_map[old] for old in df.columns if old in col_map})

    # Keep only required columns
    if 'Protein' in df.columns and 'external_score' in df.columns:
        return df[['Protein', 'external_score']]
    else:
        raise ValueError("Could not find both 'Protein' and 'external_score' (or their equivalents).")