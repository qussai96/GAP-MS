import re
import pandas as pd
import csv
from pathlib import Path


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

def save_gtf_subset_all_features(gtf_file, subset, output_dir, output_file_name):
    """
    Read original GTF with ALL features and save only those belonging to supported proteins.
    This preserves exons, UTRs, start/stop codons, etc.
    
    Optimized for large GTF files using line-by-line parsing.
    Handles both GTF formats:
    - With explicit protein_id attributes in CDS rows
    - Where mRNA ID itself IS the protein identifier
    """
    gtf_file = Path(gtf_file)
    output_dir = Path(output_dir)
    
    # Parse GTF line by line to build ID mappings
    id_to_protein = {}  # Maps feature ID -> protein ID
    id_to_type = {}  # Maps ID -> feature type
    
    # First pass: scan for mRNA/transcript IDs and CDS protein IDs
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            
            feature_type = parts[2].lower()
            attributes = parts[8]
            
            # Extract a stable feature identifier from either GFF3 or GTF style attrs
            feature_id = (
                extract_attribute(attributes, 'ID')
                or extract_attribute(attributes, 'transcript_id')
                # Bare-ID fallback: some tools (e.g. BRAKER) write the raw ID as
                # the entire column-9 value with no key=value pairs
                or (attributes.strip() if re.match(r'^[\w.\-]+$', attributes.strip()) else None)
            )
            
            if not feature_id:
                continue
            
            id_to_type[feature_id] = feature_type
            
            # For mRNA/transcript rows: ID itself is the protein identifier
            if feature_type in ['mrna', 'transcript']:
                id_to_protein[feature_id] = feature_id
            
            # For CDS rows: check for explicit protein_id
            if feature_type == 'cds':
                protein_id = (
                    extract_attribute(attributes, 'protein_id')
                    or extract_attribute(attributes, 'transcript_id')
                )
                parent_id = (
                    extract_attribute(attributes, 'Parent')
                    or extract_attribute(attributes, 'transcript_id')
                )
                if protein_id and parent_id:
                    id_to_protein[parent_id] = protein_id
    
    # Second pass: link genes to their mRNA children
    gene_to_mrna = {}  # Maps gene ID -> mRNA ID
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            
            feature_type = parts[2].lower()
            attributes = parts[8]
            
            if feature_type in ['mrna', 'transcript']:
                mrna_id = (
                    extract_attribute(attributes, 'ID')
                    or extract_attribute(attributes, 'transcript_id')
                    or (attributes.strip() if re.match(r'^[\w.\-]+$', attributes.strip()) else None)
                )
                gene_id = (
                    extract_attribute(attributes, 'Parent')
                    or extract_attribute(attributes, 'gene_id')
                )
                if mrna_id and gene_id:
                    if gene_id not in gene_to_mrna:
                        gene_to_mrna[gene_id] = mrna_id  # Use first mRNA found
    
    # Link genes to proteins
    for gene_id, mrna_id in gene_to_mrna.items():
        if mrna_id in id_to_protein:
            id_to_protein[gene_id] = id_to_protein[mrna_id]
    
    # Third pass: filter and write output
    with open(gtf_file, 'r') as infile, \
         open(output_dir / output_file_name, 'w') as outfile:
        
        for line in infile:
            if line.startswith('#'):
                continue
            
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            
            attributes = parts[8]
            feature_id = None
            parent_id = None
            protein_id = None
            
            # Try to get protein_id directly from attributes
            protein_id = (
                extract_attribute(attributes, 'protein_id')
                or extract_attribute(attributes, 'transcript_id')
            )
            
            # Try ID attribute
            if not protein_id:
                feature_id = (
                    extract_attribute(attributes, 'ID')
                    or extract_attribute(attributes, 'transcript_id')
                )
                if feature_id:
                    protein_id = id_to_protein.get(feature_id)
            
            # Try Parent attribute
            if not protein_id:
                parent_id = (
                    extract_attribute(attributes, 'Parent')
                    or extract_attribute(attributes, 'gene_id')
                )
                if parent_id:
                    protein_id = id_to_protein.get(parent_id)
            
            # Bare-ID fallback: attrs column is a plain identifier with no key=value pairs
            if not protein_id and re.match(r'^[\w.\-]+$', attributes.strip()):
                bare = attributes.strip()
                protein_id = id_to_protein.get(bare) or (bare if bare in subset else None)
            
            # Write if protein is in subset
            if protein_id and protein_id in subset:
                outfile.write(line)