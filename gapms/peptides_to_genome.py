#!/usr/bin/env python3
"""
Map peptides to genomic coordinates for visualization in genome browsers.
Outputs BED or GFF3 format files with peptide locations on the genome.
"""

from pathlib import Path
import pandas as pd
from Bio import SeqIO
from typing import List, Tuple, Dict

from .gtf_utils import gtf_to_df_with_genes
from .peptide_utils import mapping_file_to_df, get_gene_protein_specific_peps


def get_protein_to_cds_mapping(gtf_df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """
    Create a mapping from protein ID to sorted CDS features.
    
    Returns:
        Dict mapping protein ID to DataFrame of CDS features sorted by genomic position
    """
    cds_df = gtf_df[gtf_df['Type'] == 'CDS'].copy()
    
    # Sort CDS by genomic position for each protein
    cds_df = cds_df.sort_values(['Protein', 'Start'])
    
    # Calculate cumulative codon positions
    df_CDSs = cds_df.drop_duplicates(subset=['Seqid', 'Start', 'End', 'Protein']).copy()
    df_CDSs['cds_len'] = ((df_CDSs['End'] - df_CDSs['Start']) / 3).round().astype(int)
    
    # Forward strand (+)
    forward = df_CDSs['Strand'] == '+'
    df_CDSs.loc[forward, 'cds_start'] = (
        df_CDSs[forward].groupby('Protein')['cds_len'].cumsum() - df_CDSs[forward]['cds_len'] + 1
    )
    df_CDSs.loc[forward, 'cds_end'] = df_CDSs[forward].groupby('Protein')['cds_len'].cumsum()

    # Reverse strand (-)
    reverse = ~forward
    reverse_group = df_CDSs[reverse].groupby('Protein')['cds_len']
    cumsum_reverse = reverse_group.cumsum()
    max_len = reverse_group.transform('sum')
    df_CDSs.loc[reverse, 'cds_end'] = max_len - cumsum_reverse + df_CDSs.loc[reverse, 'cds_len']
    df_CDSs.loc[reverse, 'cds_start'] = df_CDSs.loc[reverse, 'cds_end'] - df_CDSs.loc[reverse, 'cds_len'] + 1
    
    # Group by protein
    protein_cds_dict = {}
    for protein in df_CDSs['Protein'].unique():
        protein_cds_dict[protein] = df_CDSs[df_CDSs['Protein'] == protein].copy()
    
    return protein_cds_dict


def map_peptide_to_codons(pep_start_aa: int, pep_end_aa: int) -> List[Tuple[int, int]]:
    """
    Convert peptide amino acid positions to codon positions.
    
    Args:
        pep_start_aa: Peptide start position in amino acids (1-based)
        pep_end_aa: Peptide end position in amino acids (1-based)
    
    Returns:
        List of (codon_start, codon_end) tuples (1-based)
    """
    # Each amino acid = 1 codon
    return [(pep_start_aa, pep_end_aa)]


def map_codon_to_genome(codon_start: int, codon_end: int, 
                        cds_features: pd.DataFrame, strand: str) -> List[Tuple[int, int]]:
    """
    Map codon positions to genomic coordinates.
    
    Args:
        codon_start: Start codon position (1-based, in protein coordinates)
        codon_end: End codon position (1-based, in protein coordinates)
        cds_features: DataFrame of CDS features for this protein
        strand: + or -
    
    Returns:
        List of (genomic_start, genomic_end) tuples (1-based genomic coordinates)
    """
    genomic_regions = []
    
    for _, cds in cds_features.iterrows():
        cds_start_codon = cds['cds_start']
        cds_end_codon = cds['cds_end']
        
        # Check if peptide overlaps with this CDS
        if codon_end < cds_start_codon or codon_start > cds_end_codon:
            continue
        
        # Calculate overlap
        overlap_start = max(codon_start, cds_start_codon)
        overlap_end = min(codon_end, cds_end_codon)
        
        # Convert codon positions to nucleotide offsets within this CDS
        if strand == '+':
            # Forward strand: codon position increases with genomic position
            nt_offset_start = (overlap_start - cds_start_codon) * 3
            nt_offset_end = (overlap_end - cds_start_codon + 1) * 3 - 1
            
            genomic_start = cds['Start'] + nt_offset_start
            genomic_end = cds['Start'] + nt_offset_end
        else:
            # Reverse strand: codon position decreases with genomic position
            nt_offset_start = (overlap_start - cds_start_codon) * 3
            nt_offset_end = (overlap_end - cds_start_codon + 1) * 3 - 1
            
            genomic_end = cds['End'] - nt_offset_start
            genomic_start = cds['End'] - nt_offset_end
        
        genomic_regions.append((int(genomic_start), int(genomic_end)))
    
    return genomic_regions


def classify_peptide_type(genomic_regions: List[Tuple[int, int]], 
                         is_first_aa: bool, is_last_aa: bool) -> str:
    """
    Classify peptide based on its genomic mapping.
    
    Returns:
        'spliced' if spans multiple exons, 'N-terminal' if at start,
        'C-terminal' if at end, 'internal' otherwise
    """
    if len(genomic_regions) > 1:
        return 'spliced'
    elif is_first_aa:
        return 'N-terminal'
    elif is_last_aa:
        return 'C-terminal'
    else:
        return 'internal'


def export_to_bed(peptide_mappings: List[Dict], output_file: Path):
    """
    Export peptide mappings to BED12 format for genome browser visualization.

    This function writes only a BED file. The "name" field contains only the peptide
    sequence (no protein or type). Spliced peptides are represented as multi-block
    BED entries (blockCount > 1).
    """
    # Filter out mappings with empty genomic regions before sorting
    valid_mappings = [m for m in peptide_mappings if m.get('genomic_regions')]
    
    # Sort mappings by seqid (chromosome) and then by start position
    sorted_mappings = sorted(valid_mappings, 
                            key=lambda m: (m['seqid'], min([r[0] for r in m['genomic_regions']])))
    
    with open(output_file, 'w') as bed:
        bed.write("track name=Peptides description=\"Mapped Peptides\" itemRgb=\"On\"\n")

        for mapping in sorted_mappings:
            chrom = mapping['seqid']
            strand = mapping['strand']
            peptide = mapping['peptide']
            regions = mapping['genomic_regions']

            if not regions:
                continue

            # Color coding by type (use mapping['type'] if available)
            colors = {
                'spliced': '255,0,0',      # Red
                'N-terminal': '0,255,0',   # Green
                'C-terminal': '51,51,51',   # Blue
                'internal': '153,153,153'  # Gray
            }
            color = colors.get(mapping.get('type', 'internal'), '0,0,0')

            # Overall start and end
            all_starts = [r[0] for r in regions]
            all_ends = [r[1] for r in regions]
            chrom_start = min(all_starts) - 1  # BED is 0-based
            chrom_end = max(all_ends)

            # Block information for spliced peptides
            block_count = len(regions)
            block_sizes = ','.join([str(end - start + 1) for start, end in regions])
            block_starts = ','.join([str(start - 1 - chrom_start) for start, _ in regions])

            # Name is peptide sequence (reversed for reverse strand for proper IGV display)
            name = peptide[::-1] if strand == '-' else peptide
            score = 0

            # BED fields: chrom, start, end, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
            bed.write(f"{chrom}\t{chrom_start}\t{chrom_end}\t{name}\t{score}\t{strand}\t")
            bed.write(f"{chrom_start}\t{chrom_end}\t{color}\t")
            bed.write(f"{block_count}\t{block_sizes}\t{block_starts}\n")


def export_to_gff3(peptide_mappings: List[Dict], output_file: Path):
    """Legacy: GFF3 export retained for compatibility but not used by default."""
    # Keep implementation minimal in case it's needed in future
    with open(output_file, 'w') as gff:
        gff.write("##gff-version 3\n")
        # Not used by default


def map_peptides_to_genome(gtf_file: Path, protein_fasta: Path, mapping_file: Path,
                           output_dir: Path):
    """
    Main function to map peptides to genomic coordinates.

    This function ONLY writes a BED12 file (peptides_mapped.bed) in the output_dir.
    The BED "name" field will contain only the peptide sequence.
    """
    print("Loading GTF annotations...")
    gtf_df = gtf_to_df_with_genes(gtf_file)
    
    print("Loading protein sequences...")
    seq_dict = {record.id: str(record.seq).replace('I', 'L').replace('*', '') 
                for record in SeqIO.parse(protein_fasta, "fasta")}
    
    print("Loading peptide mappings...")
    mapping_df = mapping_file_to_df(mapping_file)
    # Normalize peptides: replace I with L (isoleucine/leucine are isomeric in MS)
    mapping_df['Peptide'] = mapping_df['Peptide'].str.replace('I', 'L')
    pep_df = get_gene_protein_specific_peps(gtf_df, mapping_df)
    pep_df = pep_df.drop_duplicates(
    subset=["Peptide", "Protein", "Gene"], keep='first')
    print("Building protein-to-CDS mapping...")
    protein_cds_dict = get_protein_to_cds_mapping(gtf_df)
    
    print("Mapping peptides to genomic coordinates...")
    peptide_mappings = []
    skipped_not_found = []
    skipped_no_cds = []
    
    for _, row in pep_df.iterrows():
        peptide = row['Peptide']
        protein = row['Protein']
        gene = row['Gene']
        
        # Get protein sequence and peptide position
        prot_seq = seq_dict.get(protein)
        if not prot_seq:
            skipped_no_cds.append((peptide, protein, "protein not in seq_dict"))
            continue
        
        pep_start = prot_seq.find(peptide) + 1  # 1-based
        if pep_start == 0:
            skipped_not_found.append((peptide, protein, "peptide not found in protein"))
            continue
        
        pep_end = pep_start + len(peptide) - 1
        prot_len = len(prot_seq)
        
        # Get CDS features for this protein
        if protein not in protein_cds_dict:
            continue
        
        cds_features = protein_cds_dict[protein].copy()
        if cds_features.empty:
            continue
        
        seqid = cds_features.iloc[0]['Seqid']
        strand = cds_features.iloc[0]['Strand']
        
        # For selenoproteins: extend CDS if peptide is beyond annotated CDS
        # This handles cases where protein extends beyond stop codon
        max_cds_end = cds_features['cds_end'].max()
        if pep_end > max_cds_end:
            # Extend the last CDS feature to cover the full protein
            last_cds_idx = cds_features['cds_end'].idxmax()
            additional_codons = pep_end - max_cds_end
            additional_nt = additional_codons * 3
            
            # Extend genomic coordinates based on strand
            if strand == '+':
                cds_features.loc[last_cds_idx, 'End'] += additional_nt
            else:
                cds_features.loc[last_cds_idx, 'Start'] -= additional_nt
            
            # Update CDS end position
            cds_features.loc[last_cds_idx, 'cds_end'] = pep_end
            cds_features.loc[last_cds_idx, 'cds_len'] = int((cds_features.loc[last_cds_idx, 'End'] - 
                                                              cds_features.loc[last_cds_idx, 'Start']) / 3)
        
        # Map peptide to genomic coordinates
        genomic_regions = map_codon_to_genome(pep_start, pep_end, cds_features, strand)
        
        # Skip peptides with no genomic regions
        if not genomic_regions:
            skipped_no_cds.append((peptide, protein, "no genomic regions found"))
            continue
        
        # Classify peptide type
        is_first_aa = (pep_start == 1)
        is_last_aa = (pep_end == prot_len)
        pep_type = classify_peptide_type(genomic_regions, is_first_aa, is_last_aa)
        
        peptide_mappings.append({
            'seqid': seqid,
            'strand': strand,
            'peptide': peptide,
            'protein': protein,
            'gene': gene,
            'type': pep_type,
            'pep_start': pep_start,
            'pep_end': pep_end,
            'genomic_regions': genomic_regions
        })
    
    print(f"Mapped {len(peptide_mappings)} peptides to genomic coordinates")
    
    # Report on skipped peptides
    if skipped_not_found:
        print(f"\nSkipped (peptide not found in protein sequence): {len(skipped_not_found)}")
        # Show first few examples with 'I'
        for pep, prot, reason in skipped_not_found[:5]:
            if 'I' in pep:
                print(f"  {pep} in {prot}")
    
    if skipped_no_cds:
        print(f"Skipped (protein not in seq_dict): {len(skipped_no_cds)}")
    
    # Remove duplicate entries: same seqid, strand, and genomic regions
    # (even if from different proteins or genes)
    seen = set()
    unique_mappings = []
    for m in peptide_mappings:
        # Create a hashable key for the genomic location and peptide
        key = (m['seqid'], m['strand'], m['peptide'], 
               tuple(sorted(m['genomic_regions'])))
        if key not in seen:
            seen.add(key)
            unique_mappings.append(m)
    
    print(f"After deduplication: {len(unique_mappings)} unique peptide mappings "
          f"(removed {len(peptide_mappings) - len(unique_mappings)} duplicates)")
    peptide_mappings = unique_mappings
    
    # Export only BED (per user request)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    bed_file = output_dir / 'peptides_mapped.bed'
    print(f"Writing BED file: {bed_file}")
    export_to_bed(peptide_mappings, bed_file)
    
    # Print summary statistics
    type_counts = pd.Series([m['type'] for m in peptide_mappings]).value_counts()
    print("\nPeptide type summary:")
    for pep_type, count in type_counts.items():
        print(f"  {pep_type}: {count}")
    
    print("\nDone! Files ready for genome browser visualization.")


def main():
    """Entry point for command-line tool."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Map peptides to genomic coordinates for genome browser visualization"
    )
    parser.add_argument("-g", "--gtf", type=Path, required=True,
                       help="GTF file with gene predictions")
    parser.add_argument("-f", "--fasta", type=Path, required=True,
                       help="Protein FASTA file")
    parser.add_argument("-m", "--mapping", type=Path, required=True,
                       help="Peptide-to-protein mapping file")
    parser.add_argument("-o", "--output", type=Path, required=True,
                       help="Output directory")
    
    args = parser.parse_args()
    
    map_peptides_to_genome(
        args.gtf,
        args.fasta,
        args.mapping,
        args.output
    )


if __name__ == "__main__":
    main()
