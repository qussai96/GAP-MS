import re
import pandas as pd
import csv

def extract_id(suppl_value, pattern):
    if isinstance(suppl_value, str):
        match = re.search(pattern, suppl_value)
        if match:
            return match.group(1)
    # print(f'Check row with suppl {suppl_value}')
    return None

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
        suppl = row['Suppl']
        
        # Try gene_id first (GTF/GFF)
        gene_id = extract_attribute(suppl, 'gene_id')
        if gene_id:
            return gene_id
        
        # Try GeneID (NCBI GFF)
        gene_id = extract_attribute(suppl, 'GeneID')
        if gene_id:
            return gene_id
        
        # Try gene attribute
        gene_id = extract_attribute(suppl, 'gene')
        if gene_id:
            return gene_id
        
        # Try Parent for non-gene features
        if row['Type'].lower() != 'gene':
            parent = extract_attribute(suppl, 'Parent')
            if parent:
                return parent
        
        # Try ID for gene features
        if row['Type'].lower() == 'gene':
            id_val = extract_attribute(suppl, 'ID')
            if id_val:
                return id_val
        
        # Fallback for transcript/mRNA/gene
        if row['Type'].lower() in ['mrna', 'transcript', 'gene']:
            return suppl
        
        return None

    def protein_extraction(row):
        suppl = row['Suppl']
        row_type = row['Type'].lower()
        
        if row_type in ['mrna', 'transcript', 'gene']:
            # Try ID first
            id_val = extract_attribute(suppl, 'ID')
            if id_val:
                return id_val
            
            # Try transcript_id
            transcript_id = extract_attribute(suppl, 'transcript_id')
            if transcript_id:
                return transcript_id
            
            return suppl
        else:
            # CDS features - try protein_id first
            protein_id = extract_attribute(suppl, 'protein_id')
            if protein_id:
                return protein_id
            
            # Try transcript_id
            transcript_id = extract_attribute(suppl, 'transcript_id')
            if transcript_id:
                return transcript_id
            
            # Try Parent
            parent = extract_attribute(suppl, 'Parent')
            if parent:
                return parent
            
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