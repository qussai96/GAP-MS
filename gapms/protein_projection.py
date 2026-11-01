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
from .peptide_utils import (
    mapping_file_to_df,
    get_gene_protein_specific_peps,
    resolve_peptide_starts,
)


def _frame_value(value) -> int:
    """Return a valid GTF/GFF CDS phase, defaulting missing values to zero."""
    try:
        phase = int(value)
    except (TypeError, ValueError):
        return 0
    return phase if phase in {0, 1, 2} else 0


def _coding_nt_positions(cds_features: pd.DataFrame) -> List[int]:
    """Return CDS genomic positions in translation order.

    The phase on the first CDS trims an incomplete leading codon. Later phases
    describe codons split across adjacent CDS records; concatenating the records
    preserves those split codons.
    """
    if cds_features.empty:
        return []

    strand = cds_features.iloc[0]['Strand']
    ordered = cds_features.sort_values('Start', ascending=strand != '-')
    positions = []
    for index, (_, cds) in enumerate(ordered.iterrows()):
        start, end = int(cds['Start']), int(cds['End'])
        feature_positions = (
            list(range(start, end + 1))
            if strand == '+'
            else list(range(end, start - 1, -1))
        )
        if index == 0:
            feature_positions = feature_positions[_frame_value(cds.get('Frame')):]
        positions.extend(feature_positions)
    return positions


def _annotate_cds_coordinates_reference(cds_df: pd.DataFrame) -> pd.DataFrame:
    """Add translation-space nucleotide/amino-acid bounds to CDS rows."""
    annotated_groups = []
    for _, group in cds_df.groupby('Protein', sort=False, dropna=False):
        if group.empty:
            continue
        strand = group.iloc[0]['Strand']
        ordered = group.sort_values('Start', ascending=strand != '-').copy()
        cumulative_nt = 0
        nt_starts = []
        nt_ends = []
        for index, (_, cds) in enumerate(ordered.iterrows()):
            nt_length = int(cds['End']) - int(cds['Start']) + 1
            if index == 0:
                nt_length -= _frame_value(cds.get('Frame'))
            nt_length = max(nt_length, 0)
            nt_starts.append(cumulative_nt + 1)
            cumulative_nt += nt_length
            nt_ends.append(cumulative_nt)

        ordered['cds_nt_start'] = nt_starts
        ordered['cds_nt_end'] = nt_ends
        ordered['cds_start'] = ((ordered['cds_nt_start'] - 1) // 3) + 1
        ordered['cds_end'] = (ordered['cds_nt_end'] + 2) // 3
        ordered['cds_len'] = ordered['cds_nt_end'] - ordered['cds_nt_start'] + 1
        annotated_groups.append(ordered)

    if not annotated_groups:
        return cds_df.assign(
            cds_nt_start=pd.Series(dtype='Int64'),
            cds_nt_end=pd.Series(dtype='Int64'),
            cds_start=pd.Series(dtype='Int64'),
            cds_end=pd.Series(dtype='Int64'),
            cds_len=pd.Series(dtype='Int64'),
        )
    return pd.concat(annotated_groups).sort_index()


def annotate_cds_coordinates(cds_df: pd.DataFrame) -> pd.DataFrame:
    """Vectorized CDS offsets, preserving the original grouping/sort semantics."""
    import numpy as np
    if (cds_df.empty or not cds_df.index.is_unique
            or not pd.api.types.is_numeric_dtype(cds_df['Start'])
            or cds_df.duplicated(['Protein', 'Start']).any()):
        # Retain pandas' original per-group tie ordering for unusual inputs.
        return _annotate_cds_coordinates_reference(cds_df)
    codes, _ = pd.factorize(cds_df['Protein'], sort=False, use_na_sentinel=False)
    work = cds_df.copy()
    work['_cds_group'] = codes
    first_strands = work.drop_duplicates('_cds_group').set_index('_cds_group')['Strand']
    work['_cds_order'] = np.where(work['_cds_group'].map(first_strands).eq('-'),
                                  -work['Start'], work['Start'])
    work = work.sort_values(['_cds_group', '_cds_order'])
    lengths = work['End'].astype('int64') - work['Start'].astype('int64') + 1
    first = ~work['_cds_group'].duplicated()
    lengths.loc[first] -= work.loc[first, 'Frame'].map(_frame_value) if 'Frame' in work else 0
    lengths = lengths.clip(lower=0)
    ends = lengths.groupby(work['_cds_group'], sort=False).cumsum()
    work['cds_nt_start'] = ends - lengths + 1
    work['cds_nt_end'] = ends
    work['cds_start'] = (work['cds_nt_start'] - 1) // 3 + 1
    work['cds_end'] = (ends + 2) // 3
    work['cds_len'] = lengths
    return work.drop(columns=['_cds_group', '_cds_order']).sort_index()


def get_protein_to_cds_mapping(
    gtf_df: pd.DataFrame, protein_ids=None
) -> Dict[str, pd.DataFrame]:
    """
    Create a mapping from protein ID to sorted CDS features.
    
    Returns:
        Dict mapping protein ID to DataFrame of CDS features sorted by genomic position
    """
    cds_df = gtf_df[gtf_df['Type'] == 'CDS'].copy()
    if protein_ids is not None:
        protein_ids = set(protein_ids)
        cds_df = cds_df[cds_df['Protein'].isin(protein_ids)]
    
    df_CDSs = cds_df.drop_duplicates(subset=['Seqid', 'Start', 'End', 'Protein']).copy()
    df_CDSs = annotate_cds_coordinates(df_CDSs)
    
    # A boolean scan for every protein is quadratic on large annotations.
    # Group once so construction remains linear in the number of CDS rows.
    return {
        protein: group.copy()
        for protein, group in df_CDSs.groupby('Protein', sort=False)
    }


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
    if codon_start < 1 or codon_end < codon_start:
        return []

    positions = _coding_nt_positions(cds_features)
    nt_start = (codon_start - 1) * 3
    nt_end = codon_end * 3
    if nt_end > len(positions):
        return []

    selected = positions[nt_start:nt_end]
    if not selected:
        return []

    genomic_regions = []
    region = [selected[0]]
    for position in selected[1:]:
        if abs(position - region[-1]) == 1:
            region.append(position)
        else:
            genomic_regions.append((min(region), max(region)))
            region = [position]
    genomic_regions.append((min(region), max(region)))
    return sorted(genomic_regions)


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
    sorted_mappings = sorted(
        valid_mappings,
        key=lambda m: (
            str(m['seqid']),
            min(int(region[0]) for region in m['genomic_regions']),
        ),
    )
    
    with open(output_file, 'w') as bed:
        bed.write("track name=Peptides description=\"Mapped Peptides\" itemRgb=\"On\"\n")

        for mapping in sorted_mappings:
            chrom = str(mapping['seqid'])
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

            # Peptide sequences are always written N-to-C; genomic strand does
            # not reverse an amino-acid sequence.
            name = peptide
            score = 0

            # BED fields: chrom, start, end, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
            bed.write(f"{chrom}\t{chrom_start}\t{chrom_end}\t{name}\t{score}\t{strand}\t")
            bed.write(f"{chrom_start}\t{chrom_end}\t{color}\t")
            bed.write(f"{block_count}\t{block_sizes}\t{block_starts}\t")
            bed.write(
                f"protein={mapping['protein']}\tpep_start={mapping['pep_start']}\t"
                f"pep_end={mapping['pep_end']}\tprot_len={mapping['prot_len']}\t"
                f"previous_aa={mapping['previous_aa']}\n"
            )


def map_peptides_to_genome(gtf_file: Path, protein_fasta: Path, mapping_file: Path,
                           output_dir: Path):
    """
    Main function to map peptides to genomic coordinates.

    This function ONLY writes a BED12 file (peptides_mapped.bed) in the output_dir.
    The BED "name" field will contain only the peptide sequence.
    """
    print("Loading protein sequences...")
    seq_dict = {record.id: str(record.seq).replace('I', 'L').replace('*', '') 
                for record in SeqIO.parse(protein_fasta, "fasta")}

    print("Loading GTF annotations...")
    gtf_df = gtf_to_df_with_genes(gtf_file, valid_protein_ids=seq_dict)
    
    print("Loading peptide mappings...")
    mapping_df = mapping_file_to_df(mapping_file)
    # Normalize peptides: replace I with L (isoleucine/leucine are isomeric in MS)
    mapping_df['Peptide'] = mapping_df['Peptide'].str.replace('I', 'L')
    pep_df = get_gene_protein_specific_peps(gtf_df, mapping_df)
    dedupe_columns = ["Peptide", "Protein", "Gene"]
    if 'location' in pep_df.columns:
        dedupe_columns.append('location')
    pep_df = pep_df.drop_duplicates(subset=dedupe_columns, keep='first')
    print("Building protein-to-CDS mapping...")
    protein_cds_dict = get_protein_to_cds_mapping(
        gtf_df, protein_ids=pep_df['Protein'].dropna().unique()
    )
    
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
        
        pep_starts = resolve_peptide_starts(
            prot_seq,
            peptide,
            row.get('location'),
        )
        if not pep_starts:
            skipped_not_found.append((peptide, protein, "peptide not found in protein"))
            continue

        # Get CDS features for this protein
        if protein not in protein_cds_dict:
            continue

        cds_features = protein_cds_dict[protein].copy()
        if cds_features.empty:
            continue

        seqid = cds_features.iloc[0]['Seqid']
        strand = cds_features.iloc[0]['Strand']

        prot_len = len(prot_seq)
        for pep_start in pep_starts:
            pep_end = pep_start + len(peptide) - 1

            # Map peptide to genomic coordinates
            genomic_regions = map_codon_to_genome(
                pep_start, pep_end, cds_features, strand
            )

            # Skip peptides with no genomic regions
            if not genomic_regions:
                skipped_no_cds.append((peptide, protein, "no genomic regions found"))
                continue

            # Classify peptide type
            is_first_aa = (pep_start == 1)
            is_last_aa = (pep_end == prot_len)
            pep_type = classify_peptide_type(
                genomic_regions, is_first_aa, is_last_aa
            )

            peptide_mappings.append({
                'seqid': seqid,
                'strand': strand,
                'peptide': peptide,
                'protein': protein,
                'gene': gene,
                'type': pep_type,
                'pep_start': pep_start,
                'pep_end': pep_end,
                'prot_len': prot_len,
                'previous_aa': prot_seq[pep_start - 2] if pep_start > 1 else '',
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
    
    # Remove only true duplicate mappings. Keep the same genomic peptide once
    # per protein because downstream event calls require protein-specific
    # evidence and otherwise one isoform could incorrectly consume another's
    # terminal or junction peptide.
    seen = set()
    unique_mappings = []
    for m in peptide_mappings:
        # Create a hashable key for the genomic location and peptide
        key = (
            m['seqid'], m['strand'], m['peptide'], m['protein'],
            m['pep_start'], m['pep_end'], tuple(sorted(m['genomic_regions'])),
        )
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
