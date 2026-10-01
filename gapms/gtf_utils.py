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


def _build_identifier_lookup(valid_ids):
    """Build exact and unambiguous version-insensitive FASTA ID lookups."""
    exact = {str(value) for value in valid_ids if value is not None}
    aliases = {}
    ambiguous = set()

    for value in exact:
        unversioned = re.sub(r'\.\d+$', '', value)
        if unversioned in aliases and aliases[unversioned] != value:
            ambiguous.add(unversioned)
        else:
            aliases[unversioned] = value

    for value in ambiguous:
        aliases.pop(value, None)
    return exact, aliases


def _resolve_identifier(value, exact_ids, aliases):
    """Return the corresponding FASTA ID for an annotation identifier."""
    if value is None or pd.isna(value):
        return None
    value = str(value)
    if value in exact_ids:
        return value
    return aliases.get(re.sub(r'\.\d+$', '', value))


def _select_matching_identifier(candidates, valid_ids, fallback):
    """Select the first candidate present in the FASTA, else use fallback."""
    if valid_ids is None:
        return fallback

    exact_ids, aliases = _build_identifier_lookup(valid_ids)
    selected = pd.Series(index=fallback.index, dtype='object')
    for candidate in candidates:
        matched = candidate.map(
            lambda value: _resolve_identifier(value, exact_ids, aliases)
        )
        selected = selected.combine_first(matched)
    return selected.combine_first(fallback)


def gtf_to_df_with_genes(gtf_file, valid_protein_ids=None):
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
        # Sequence identifiers are labels, even when every value in one input
        # happens to be numeric. Keeping them as strings prevents mixed-type
        # contig collections from failing during deterministic BED sorting.
        dtype={'Seqid': str, 'Type': str, 'Suppl': str},
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

    id_attr = _str_extract(s, r'ID=([^;,\s]+)')
    transcript_attr = _str_extract(s, r'transcript_id\s*[= ]?"?([^";,\s]+)"?')
    protein_attr = _str_extract(s, r'protein_id\s*[= ]?"?([^";,\s]+)"?')
    parent_attr = _str_extract(s, r'Parent=([^;,\s]+)')

    # CDS identifiers are chosen from the namespace actually present in the
    # FASTA. This handles annotations that provide both protein_id and
    # transcript_id while gffread names translated records by transcript_id.
    cds_mask = t == 'cds'
    cds_fallback = (
        protein_attr.where(cds_mask)
        .combine_first(transcript_attr.where(cds_mask))
        .combine_first(parent_attr.where(cds_mask))
    )
    prot_cds = _select_matching_identifier(
        [
            protein_attr.where(cds_mask),
            transcript_attr.where(cds_mask),
            parent_attr.where(cds_mask),
        ],
        valid_protein_ids,
        cds_fallback,
    )

    # Link transcript rows to the identifier selected for their CDS children.
    # This preserves transcript-level scores when the FASTA uses protein IDs.
    transcript_to_protein = {}
    if valid_protein_ids is not None:
        links = pd.DataFrame({
            'transcript': transcript_attr.where(cds_mask).combine_first(parent_attr.where(cds_mask)),
            'protein': prot_cds,
        }).dropna().drop_duplicates()
        if not links.empty:
            unique_links = links.groupby('transcript')['protein'].nunique()
            unique_transcripts = unique_links[unique_links == 1].index
            transcript_to_protein = (
                links[links['transcript'].isin(unique_transcripts)]
                .drop_duplicates('transcript')
                .set_index('transcript')['protein']
                .to_dict()
            )

    # Vectorized protein extraction for transcript/mRNA/gene rows.
    meta_fallback = (
        id_attr.where(meta_mask)
        .combine_first(transcript_attr.where(meta_mask))
        .combine_first(s.where(meta_mask))
    )
    linked_meta = transcript_attr.where(meta_mask).map(transcript_to_protein)
    prot_meta = _select_matching_identifier(
        [
            linked_meta,
            id_attr.where(meta_mask),
            transcript_attr.where(meta_mask),
            protein_attr.where(meta_mask),
        ],
        valid_protein_ids,
        meta_fallback,
    )
    df['Protein'] = prot_meta.combine_first(prot_cds)

    # Normalize Type values to what downstream code expects
    type_map = {'gene': 'gene', 'mrna': 'mRNA', 'transcript': 'transcript', 'cds': 'CDS'}
    mapped_type = df['Type'].map(type_map)
    df['Type'] = mapped_type.where(mapped_type.notna(), df['Type'])

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
    subset = {str(protein_id) for protein_id in subset}
    exact_ids, aliases = _build_identifier_lookup(subset)

    def matching_id(*candidates):
        for candidate in candidates:
            resolved = _resolve_identifier(candidate, exact_ids, aliases)
            if resolved is not None:
                return resolved
        return None
    
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
                id_to_protein[feature_id] = matching_id(feature_id) or feature_id
            
            # For CDS rows: check for explicit protein_id
            if feature_type == 'cds':
                protein_id = matching_id(
                    extract_attribute(attributes, 'protein_id'),
                    extract_attribute(attributes, 'transcript_id'),
                    extract_attribute(attributes, 'Parent'),
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
            
            # Prefer whichever annotation identifier matches the selected
            # FASTA namespace (protein_id, transcript_id, ID, or Parent).
            protein_id = matching_id(
                extract_attribute(attributes, 'protein_id'),
                extract_attribute(attributes, 'transcript_id'),
                extract_attribute(attributes, 'ID'),
                extract_attribute(attributes, 'Parent'),
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
