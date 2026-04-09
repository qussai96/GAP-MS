#!/usr/bin/env python3
"""
Generate a TSV report identifying genes with peptide evidence at specific
locations that differ from reference annotations.

This script focuses on three categories:
- peptide_support_different_start: Genes with different start positions AND N-terminal peptide evidence
- peptide_support_different_stop: Genes with different stop positions AND C-terminal peptide evidence
- peptide_support_different_splice: Genes with different splice sites AND peptides spanning the DIFFERENT junctions

Note: For splice junctions, we verify that peptides actually span the specific junctions that 
differ from the reference, not just any splice junctions in the protein.
"""

import sys
import pandas as pd
from pathlib import Path
from collections import defaultdict
import re
from Bio import SeqIO


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


def parse_gtf_to_cds_dict(gtf_file):
    """
    Parse GTF file and return dictionary of CDS regions keyed by transcript/protein ID.
    Returns: {transcript_id: {'seqid': chr, 'strand': strand, 'cds_regions': [(start, end), ...], 'gene_id': id}}
    """
    cds_dict = defaultdict(lambda: {'seqid': None, 'strand': None, 'cds_regions': [], 'gene_id': None})
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            seqid = fields[0]
            feature_type = fields[2].lower()
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]
            
            # Only process CDS features
            if feature_type != 'cds':
                continue
            
            # Extract transcript ID (Parent)
            transcript_id = extract_attribute(attributes, 'Parent')
            if not transcript_id:
                transcript_id = extract_attribute(attributes, 'transcript_id')
            if not transcript_id:
                continue
            
            # Extract gene ID
            gene_id = extract_attribute(attributes, 'gene_id')
            if not gene_id:
                gene_id = extract_attribute(attributes, 'GeneID')
            
            cds_dict[transcript_id]['seqid'] = seqid
            cds_dict[transcript_id]['strand'] = strand
            cds_dict[transcript_id]['gene_id'] = gene_id
            cds_dict[transcript_id]['cds_regions'].append((start, end))
    
    # Sort CDS regions for each transcript
    for transcript_id in cds_dict:
        cds_dict[transcript_id]['cds_regions'].sort()
    
    return cds_dict


def get_first_and_last_cds(cds_list, strand):
    """
    Get the first and last CDS boundaries based on strand direction.
    
    For positive strand: 
        - First CDS start = leftmost start (translation start)
        - Last CDS end = rightmost end (translation end)
    For negative strand:
        - First CDS start = rightmost end (translation starts from here)
        - Last CDS end = leftmost start (translation ends here)
    
    Returns: (first_pos, last_pos) where first_pos is where translation starts,
             last_pos is where translation ends
    """
    if not cds_list:
        return None, None
    
    sorted_cds = sorted(cds_list, key=lambda x: x[0])
    
    if strand == '+':
        first_pos = sorted_cds[0][0]   # leftmost start
        last_pos = sorted_cds[-1][1]   # rightmost end
    else:
        first_pos = sorted_cds[-1][1]  # rightmost end
        last_pos = sorted_cds[0][0]    # leftmost start
    
    return first_pos, last_pos


def classify_difference(pred_info, ref_info):
    """
    Classify the type of difference between predicted and reference CDS coordinates.
    
    Returns: classification string describing the difference
    """
    if not pred_info or not ref_info:
        return "missing_coordinates"
    
    pred_cds = pred_info.get('cds_regions', [])
    ref_cds = ref_info.get('cds_regions', [])
    strand = pred_info.get('strand') or ref_info.get('strand')
    
    if not pred_cds or not ref_cds:
        return "missing_cds_regions"
    
    # Check if same seqid and strand
    if pred_info.get('seqid') != ref_info.get('seqid') or pred_info.get('strand') != ref_info.get('strand'):
        return "different_chromosome_or_strand"
    
    pred_first, pred_last = get_first_and_last_cds(pred_cds, strand)
    ref_first, ref_last = get_first_and_last_cds(ref_cds, strand)
    
    # Count exons
    pred_exons = len(pred_cds)
    ref_exons = len(ref_cds)
    
    same_first = (pred_first == ref_first)
    same_last = (pred_last == ref_last)
    
    # Classify based on CDS structure
    if same_first and same_last:
        return "same_start_same_stop" if pred_exons == ref_exons else "same_start_same_stop_different_exons"
    elif same_first and not same_last:
        return "same_start_different_stop"
    elif not same_first and same_last:
        return "different_start_same_stop"
    elif pred_first < ref_first if strand == '+' else pred_first > ref_first:
        # Predicted is upstream of reference
        return "upstream_extension"
    elif pred_first > ref_first if strand == '+' else pred_first < ref_first:
        # Predicted is downstream of reference
        return "downstream_truncation"
    else:
        if pred_exons > ref_exons:
            return "more_exons"
        elif pred_exons < ref_exons:
            return "fewer_exons"
        else:
            return "exon_structure_difference"


def get_splice_junctions(cds_regions):
    """
    Extract splice junctions from CDS regions.
    A junction is defined by the end of one exon and the start of the next exon.
    
    Returns: set of tuples (exon_end, next_exon_start)
    """
    if not cds_regions or len(cds_regions) < 2:
        return set()
    
    sorted_cds = sorted(cds_regions, key=lambda x: x[0])
    junctions = set()
    
    for i in range(len(sorted_cds) - 1):
        exon_end = sorted_cds[i][1]
        next_exon_start = sorted_cds[i + 1][0]
        junctions.add((exon_end, next_exon_start))
    
    return junctions


def load_peptides_bed(bed_file):
    """
    Load peptides_mapped.bed file and return peptide locations.
    
    Returns: dict mapping peptide_name to list of genomic locations
        {peptide_name: [{'seqid': chr, 'start': pos, 'end': pos, 'strand': strand, 'protein': protein_id}, ...]}
    """
    peptide_locations = defaultdict(list)
    
    if not Path(bed_file).exists():
        print(f"Warning: BED file not found: {bed_file}")
        return peptide_locations
    
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 6:
                continue
            
            seqid = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3]  # peptide sequence (may be reversed for minus strand)
            score = fields[4]
            strand = fields[5]
            
            # Extract protein ID from additional fields if available
            protein_id = None
            if len(fields) > 6:
                # Try to extract protein ID from attributes
                for attr in fields[6:]:
                    if 'protein' in attr.lower() or 'parent' in attr.lower():
                        protein_id = attr.split('=')[-1].split(';')[0]
                        break
            
            peptide_locations[name].append({
                'seqid': seqid,
                'start': start,
                'end': end,
                'strand': strand,
                'protein': protein_id
            })
    
    return peptide_locations


def peptide_spans_junction(peptide_loc, junction):
    """
    Check if a peptide spans a specific splice junction.
    
    Args:
        peptide_loc: dict with 'seqid', 'start', 'end', 'strand'
        junction: tuple (exon_end, next_exon_start)
    
    Returns: True if peptide spans the junction
    """
    exon_end, next_exon_start = junction
    pep_start = peptide_loc['start']
    pep_end = peptide_loc['end']
    
    # A peptide spans a junction if it overlaps both the end of one exon
    # and the start of the next exon (i.e., it crosses the intron boundary)
    # This means: pep_start <= exon_end and pep_end >= next_exon_start
    return pep_start <= exon_end and pep_end >= next_exon_start


def get_peptides_supporting_different_junctions(protein_id, pred_info, ref_info, peptide_bed_dict):
    """
    Find peptides that span splice junctions that differ between predicted and reference.
    
    Returns: number of peptides supporting the different junctions
    """
    if not pred_info or not ref_info:
        return 0
    
    pred_cds = pred_info.get('cds_regions', [])
    ref_cds = ref_info.get('cds_regions', [])
    seqid = pred_info.get('seqid')
    
    if not pred_cds or not ref_cds or len(pred_cds) < 2:
        return 0
    
    # Get junctions from both annotations
    pred_junctions = get_splice_junctions(pred_cds)
    ref_junctions = get_splice_junctions(ref_cds)
    
    # Find junctions that are different (in predicted but not in reference)
    different_junctions = pred_junctions - ref_junctions
    
    if not different_junctions:
        return 0
    
    # Find peptides for this protein
    supporting_peptides = set()
    
    # Search through all peptides in the BED file
    for peptide_name, locations in peptide_bed_dict.items():
        for loc in locations:
            # Check if this peptide location is on the same chromosome and matches protein
            if loc['seqid'] == seqid:
                # Check if peptide spans any of the different junctions
                for junction in different_junctions:
                    if peptide_spans_junction(loc, junction):
                        supporting_peptides.add(peptide_name)
                        break
    
    return len(supporting_peptides)


def generate_annotation_report(output_dir, supported_gtf, reference_gtf, 
                               tmap_file=None, all_proteins_scores=None, peptides_bed=None,
                               protein_fasta=None):
    """
    Generate comprehensive annotation comparison report.
    
    Classifies transcripts by comparing predicted GTF with reference GTF using gffcompare statistics.
    Novel genes include those with class codes: '-' (unknown), 'u' (intergenic), 'x' (antisense).
    
    Args:
        output_dir: Path to output directory
        supported_gtf: Path to supported proteins GTF file
        reference_gtf: Path to reference GTF file
        tmap_file: Path to gffcompare tmap file (auto-detected if None)
        all_proteins_scores: Path to all_proteins_scores.tsv (auto-detected if None)
        peptides_bed: Path to peptides_mapped.bed file (auto-detected if None)
        protein_fasta: Path to protein FASTA file (required for sequence export)
    """
    output_dir = Path(output_dir)
    
    # Auto-detect file paths if not provided
    if tmap_file is None:
        tmap_files = list(output_dir.glob("*.tmap"))
        if tmap_files:
            tmap_file = tmap_files[0]
        else:
            tmap_file = output_dir / "gffcompare.supported_proteins.gtf.tmap"
    
    if all_proteins_scores is None:
        all_proteins_scores = output_dir / "all_proteins_scores.tsv"
    
    if peptides_bed is None:
        peptides_bed = output_dir / "peptides_mapped.bed"
    
    tmap_file = Path(tmap_file)
    all_proteins_scores = Path(all_proteins_scores)
    peptides_bed = Path(peptides_bed)
    supported_gtf = Path(supported_gtf)
    reference_gtf = Path(reference_gtf)
    
    print(f"Loading TMAP file: {tmap_file}")
    tmap_df = pd.read_csv(tmap_file, sep='\t')
    
    print(f"Parsing GTF files...")
    supported_cds = parse_gtf_to_cds_dict(supported_gtf)
    reference_cds = parse_gtf_to_cds_dict(reference_gtf)
    
    print(f"Loading peptides BED file: {peptides_bed}")
    peptide_bed_dict = load_peptides_bed(peptides_bed)
    print(f"Loaded {len(peptide_bed_dict)} unique peptides from BED file")
    
    print(f"Loading protein scores: {all_proteins_scores}")
    if all_proteins_scores.exists():
        scores_df = pd.read_csv(all_proteins_scores, sep='\t')
        # Create dictionary with peptide type info
        scores_dict = {}
        for _, row in scores_df.iterrows():
            protein_id = row['Protein']
            scores_dict[protein_id] = {
                'mapped_peptides': row.get('mapped_peptides', 0),
                'sequence_coverage': row.get('sequence_coverage', row.get('protein' + '_coverage', 0)),
                'N_terminal_peptides': row.get('N_terminal_peptides', 0),
                'C_terminal_peptides': row.get('C_terminal_peptides', 0),
                'splice_peptides': row.get('splice_peptides', 0)
            }
    else:
        print(f"Warning: Scores file not found at {all_proteins_scores}")
        scores_dict = {}
    
    # Initialize report data
    report_data = []
    category_counts = defaultdict(lambda: {'count': 0, 'ids': set()})
    
    # Iterate through TMAP entries
    for _, row in tmap_df.iterrows():
        qry_id = row['qry_id']  # Predicted protein ID
        ref_id = row['ref_id']  # Reference protein/gene ID
        class_code = row['class_code']
        
        # Get peptide support metrics
        if qry_id in scores_dict:
            peptides = scores_dict[qry_id].get('mapped_peptides', 0)
            coverage = scores_dict[qry_id].get('sequence_coverage', 0)
            n_term_peptides = scores_dict[qry_id].get('N_terminal_peptides', 0)
            c_term_peptides = scores_dict[qry_id].get('C_terminal_peptides', 0)
            splice_peptides = scores_dict[qry_id].get('splice_peptides', 0)
        else:
            peptides = 0
            coverage = 0
            n_term_peptides = 0
            c_term_peptides = 0
            splice_peptides = 0
        
        # Skip identical matches and novel genes - we only care about differences
        if class_code == '=':
            continue
        
        # Skip novel genes
        novel_codes = ['s', 'x', 'i', 'y', 'p', 'u']
        if class_code in novel_codes:
            continue
        
        # Get predicted and reference CDS info
        pred_info = supported_cds.get(qry_id)
        ref_info = reference_cds.get(ref_id)
        
        # Classify the difference
        classification = classify_difference(pred_info, ref_info)
        
        # Determine which peptide-evidence category this belongs to
        category = None
        
        # Categories where start position differs
        start_diff_categories = ['different_start_same_stop', 'upstream_extension', 'downstream_truncation']
        # Categories where stop position differs
        stop_diff_categories = ['same_start_different_stop']
        # Categories where splice sites differ
        splice_diff_categories = ['same_start_same_stop_different_exons', 'exon_structure_difference', 
                                  'more_exons', 'fewer_exons']
        
        # Check if gene has different start AND N-terminal peptide
        if classification in start_diff_categories and n_term_peptides > 0:
            category = 'peptide_support_different_start'
        
        # Check if gene has different stop AND C-terminal peptide
        elif classification in stop_diff_categories and c_term_peptides > 0:
            category = 'peptide_support_different_stop'
        
        # Check if gene has different splice sites AND peptides spanning the DIFFERENT junctions
        elif classification in splice_diff_categories:
            # FIXED LOGIC: Check if peptides actually span the different splice junctions
            splice_support_count = get_peptides_supporting_different_junctions(
                qry_id, pred_info, ref_info, peptide_bed_dict
            )
            if splice_support_count > 0:
                category = 'peptide_support_different_splice'
                # Update splice_peptides to reflect peptides supporting DIFFERENT junctions
                splice_peptides = splice_support_count
        
        # Only include if it matches one of our 3 categories
        if category:
            # Update category counts
            if category not in category_counts:
                category_counts[category] = {'count': 0, 'ids': set()}
            
            category_counts[category]['count'] += 1
            category_counts[category]['ids'].add(qry_id)
            
            # Add to report
            report_data.append({
                'predicted_id': qry_id,
                'reference_id': ref_id,
                'structural_difference': classification,
                'peptide_evidence_category': category,
                'N_terminal_peptides': n_term_peptides,
                'C_terminal_peptides': c_term_peptides,
                'splice_peptides': splice_peptides,
                'mapped_peptides': peptides,
                'sequence_coverage': coverage
            })
    
    # Generate TSV report
    report_df = pd.DataFrame(report_data)
    
    print(f"\n{'='*80}")
    print(f"ANNOTATION COMPARISON REPORT - PEPTIDE EVIDENCE")
    print(f"{'='*80}")
    print(f"Total genes with peptide evidence at difference sites: {len(report_df)}")
    
    # Write summary report
    summary_report = []
    
    for category in sorted(category_counts.keys(), key=lambda x: category_counts[x]['count'], reverse=True):
        count = category_counts[category]['count']
        gene_ids = category_counts[category]['ids']
        
        ids_str = '{' + ','.join(sorted(list(gene_ids)[:10])) + (',...}' if len(gene_ids) > 10 else '}')
        
        print(f"{category}:")
        print(f"  Total: {count}")
        print(f"  Examples: {ids_str}")
        print()
        
        summary_report.append({
            'Category': category,
            'Count': count,
            'Protein_IDs': '{' + ','.join(sorted(list(gene_ids))) + '}'
        })
    
    # Write output files
    compare_dir = output_dir / "Compare_to_Reference"
    compare_dir.mkdir(parents=True, exist_ok=True)

    report_file = compare_dir / "annotation_comparison_report.tsv"
    summary_file = compare_dir / "annotation_comparison_summary.tsv"
    
    report_df.to_csv(report_file, sep='\t', index=False)
    print(f"\nDetailed report written to: {report_file}")
    
    summary_df = pd.DataFrame(summary_report)
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"Summary report written to: {summary_file}")
    
    # Extract GTF entries by category
    print(f"\nExtracting GTF entries by category...")
    gtf_dir = compare_dir
    different_dir = compare_dir / "Different"
    different_dir.mkdir(parents=True, exist_ok=True)
    
    # Load the supported GTF file
    supported_gtf_path = Path(supported_gtf)
    gtf_data = []
    with open(supported_gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 9:
                gtf_data.append({
                    'line': line.strip(),
                    'seqid': fields[0],
                    'source': fields[1],
                    'feature': fields[2],
                    'start': int(fields[3]),
                    'end': int(fields[4]),
                    'score': fields[5],
                    'strand': fields[6],
                    'frame': fields[7],
                    'attributes': fields[8]
                })
    
    # Load protein sequences from the provided protein FASTA
    protein_sequences = {}
    if not protein_fasta:
        raise ValueError("protein_fasta is required to export category .faa files")

    fasta_file = Path(protein_fasta)
    if fasta_file.exists():
        print(f"Loading protein sequences: {fasta_file}")
        with open(fasta_file) as fasta_handle:
            for record in SeqIO.parse(fasta_handle, "fasta"):
                protein_id = record.id.split()[0]  # Use ID before any whitespace
                protein_sequences[protein_id] = str(record.seq)
        print(f"Loaded {len(protein_sequences)} protein sequences")
    else:
        raise FileNotFoundError(f"Protein FASTA not found: {fasta_file}")
    
    # For each category, extract GTF entries and protein sequences
    for category in sorted(category_counts.keys(), key=lambda x: category_counts[x]['count'], reverse=True):
        protein_ids = category_counts[category]['ids']
        
        # Extract GTF lines for proteins in this category
        category_gtf_lines = []
        for entry in gtf_data:
            # Extract protein ID from attributes
            attrs = entry['attributes']
            protein_id = None
            
            # Try to extract protein ID from Parent attribute
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('Parent='):
                    protein_id = attr.replace('Parent=', '')
                    break
            
            if protein_id in protein_ids:
                category_gtf_lines.append(entry['line'])
        
        # Write GTF file for this category
        if category_gtf_lines:
            # Sanitize filename
            safe_category = category.replace('/', '_').replace(' ', '_')
            category_dir = different_dir if safe_category.startswith("peptide_support_different_") else gtf_dir
            gtf_file = category_dir / f"{safe_category}.gtf"
            faa_file = category_dir / f"{safe_category}.faa"
            
            with open(gtf_file, 'w') as f:
                for line in category_gtf_lines:
                    f.write(line + '\n')
            
            # Write corresponding protein sequences to .faa file
            faa_count = 0
            with open(faa_file, 'w') as f:
                for protein_id in sorted(protein_ids):
                    if protein_id in protein_sequences:
                        seq = protein_sequences[protein_id]
                        f.write(f">{protein_id}\n")
                        # Write sequence in 60-character lines (standard FASTA format)
                        for i in range(0, len(seq), 60):
                            f.write(seq[i:i+60] + '\n')
                        faa_count += 1
            
            print(f" {safe_category}.gtf ({len(category_gtf_lines)} entries) and {safe_category}.faa ({faa_count} sequences)")
    
    return report_df, summary_df


def main():
    if len(sys.argv) < 4:
        print("Usage: python compare_supp_ref.py <output_dir> <supported_gtf> <reference_gtf> [tmap_file] [scores_file] [peptides_bed] [protein_fasta]")
        print("\nArgs:")
        print("  output_dir:     Directory containing gffcompare output files")
        print("  supported_gtf:  Path to predicted GTF file (annotations)")
        print("  reference_gtf:  Path to reference GTF file")
        print("  tmap_file:      (Optional) Path to .tmap file")
        print("  scores_file:    (Optional) Path to all_proteins_scores.tsv")
        print("  peptides_bed:   (Optional) Path to peptides_mapped.bed")
        print("  protein_fasta:  (Optional) Path to protein FASTA (fallback if supported_proteins.fa is missing)")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    supported_gtf = sys.argv[2]
    reference_gtf = sys.argv[3]
    tmap_file = sys.argv[4] if len(sys.argv) > 4 else None
    scores_file = sys.argv[5] if len(sys.argv) > 5 else None
    peptides_bed = sys.argv[6] if len(sys.argv) > 6 else None
    protein_fasta = sys.argv[7] if len(sys.argv) > 7 else None
    
    generate_annotation_report(output_dir, supported_gtf, reference_gtf, tmap_file, scores_file, peptides_bed, protein_fasta)


if __name__ == '__main__':
    main()
