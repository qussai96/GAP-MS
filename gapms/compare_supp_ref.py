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
from bisect import bisect_right
import re
from Bio import SeqIO


def _pick_tmap_file(candidates):
    existing = [Path(path) for path in candidates if Path(path).exists()]
    if not existing:
        return None

    for path in existing:
        try:
            with path.open('r') as handle:
                next(handle)
                next(handle)
            return path
        except StopIteration:
            continue

    return existing[0]


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
    cds_dict = defaultdict(
        lambda: {
            'seqid': None, 'strand': None, 'cds_regions': [],
            'gene_id': None, 'protein_id': None,
        }
    )
    
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
            
            parent_id = extract_attribute(attributes, 'Parent')
            transcript_id = extract_attribute(attributes, 'transcript_id')
            protein_id = extract_attribute(attributes, 'protein_id')
            aliases = list(dict.fromkeys(
                value for value in (parent_id, transcript_id, protein_id) if value
            ))
            if not aliases:
                continue
            
            # Extract gene ID
            gene_id = extract_attribute(attributes, 'gene_id')
            if not gene_id:
                gene_id = extract_attribute(attributes, 'GeneID')
            if not gene_id:
                gene_id = extract_attribute(attributes, 'gene')
            
            preferred_protein_id = protein_id or transcript_id or parent_id
            for alias in aliases:
                cds_dict[alias]['seqid'] = seqid
                cds_dict[alias]['strand'] = strand
                cds_dict[alias]['gene_id'] = gene_id
                cds_dict[alias]['protein_id'] = preferred_protein_id
                cds_dict[alias]['cds_regions'].append((start, end))
    
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
    pred_junctions = get_splice_junctions(pred_cds)
    ref_junctions = get_splice_junctions(ref_cds)
    
    # Count exons
    pred_exons = len(pred_cds)
    ref_exons = len(ref_cds)

    # Junction differences are the strongest evidence for a splice change.
    if pred_junctions != ref_junctions:
        return "different_splice"
    
    same_first = (pred_first == ref_first)
    same_last = (pred_last == ref_last)
    last_boundary_diff = abs(pred_last - ref_last) if pred_last is not None and ref_last is not None else None
    
    # Classify based on CDS structure
    if same_first and same_last:
        return "same_start_same_stop" if pred_exons == ref_exons else "same_start_same_stop_different_exons"
    elif same_first and not same_last and (last_boundary_diff is None or last_boundary_diff > 3):
        return "same_start_different_stop"
    elif same_first and not same_last:
        return "same_start_same_stop"
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


def has_stop_difference(pred_info, ref_info, tolerance_nt=3):
    """
    Return True when translation start is unchanged but translation stop differs.

    This check is independent of splice classification so mixed splice+stop
    changes can still be categorized as stop differences.
    """
    if not pred_info or not ref_info:
        return False

    if pred_info.get('seqid') != ref_info.get('seqid') or pred_info.get('strand') != ref_info.get('strand'):
        return False

    pred_cds = pred_info.get('cds_regions', [])
    ref_cds = ref_info.get('cds_regions', [])
    if not pred_cds or not ref_cds:
        return False

    strand = pred_info.get('strand') or ref_info.get('strand')
    pred_first, pred_last = get_first_and_last_cds(pred_cds, strand)
    ref_first, ref_last = get_first_and_last_cds(ref_cds, strand)

    if pred_first != ref_first:
        return False

    if pred_last is None or ref_last is None:
        return False

    return abs(pred_last - ref_last) > tolerance_nt


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


def _donor_acceptor(junction, strand):
    """Return transcript-direction donor and acceptor CDS boundaries."""
    left_end, right_start = junction
    if strand == '-':
        return right_start, left_end
    return left_end, right_start


def build_reference_event_index(reference_cds):
    """Collect coding boundaries from every reference isoform."""
    index = defaultdict(
        lambda: {
            'starts': set(), 'stops': set(), 'donors': set(),
            'acceptors': set(), 'cds_regions': [], 'cds_starts': [],
        }
    )
    for info in reference_cds.values():
        seqid = info.get('seqid')
        strand = info.get('strand')
        cds_regions = info.get('cds_regions', [])
        if not seqid or strand not in {'+', '-'} or not cds_regions:
            continue
        key = (seqid, strand)
        first, last = get_first_and_last_cds(cds_regions, strand)
        index[key]['starts'].add(first)
        index[key]['stops'].add(last)
        index[key]['cds_regions'].extend(cds_regions)
        for junction in get_splice_junctions(cds_regions):
            donor, acceptor = _donor_acceptor(junction, strand)
            index[key]['donors'].add(donor)
            index[key]['acceptors'].add(acceptor)
    for values in index.values():
        merged_regions = []
        for start, end in sorted(values['cds_regions']):
            if merged_regions and start <= merged_regions[-1][1] + 1:
                merged_regions[-1] = (merged_regions[-1][0], max(merged_regions[-1][1], end))
            else:
                merged_regions.append((start, end))
        values['cds_regions'] = merged_regions
        values['cds_starts'] = [start for start, _ in merged_regions]
    return index


def _has_same_strand_reference_cds_overlap(pred_info, reference_event_index):
    if not pred_info:
        return False
    key = (pred_info.get('seqid'), pred_info.get('strand'))
    reference_values = reference_event_index.get(key, {})
    reference_regions = reference_values.get('cds_regions', [])
    reference_starts = reference_values.get('cds_starts', [])
    for pred_start, pred_end in pred_info.get('cds_regions', []):
        index = bisect_right(reference_starts, pred_end) - 1
        if index >= 0 and reference_regions[index][1] >= pred_start:
            return True
    return False


def _location_matches_protein_and_model(location, protein_id, pred_info):
    return (
        location.get('protein') == protein_id
        and location.get('seqid') == pred_info.get('seqid')
        and location.get('strand') == pred_info.get('strand')
    )


def build_protein_peptide_index(peptide_bed_dict):
    """Index peptide names by mapped protein without losing global ambiguity."""
    index = defaultdict(set)
    for peptide_name, locations in peptide_bed_dict.items():
        for location in locations:
            protein_id = location.get('protein')
            if protein_id:
                index[protein_id].add(peptide_name)
    return index


def get_terminal_supporting_peptides(
    protein_id, pred_info, peptide_bed_dict, terminal,
    protein_peptide_index=None,
):
    """Return peptides whose every mapping supports the same coding terminus."""
    if not pred_info or terminal not in {'start', 'stop'}:
        return set()

    strand = pred_info.get('strand')
    pred_start, pred_stop = get_first_and_last_cds(pred_info.get('cds_regions', []), strand)
    target = pred_start if terminal == 'start' else pred_stop
    supporting = set()

    candidate_peptides = (
        protein_peptide_index.get(protein_id, set())
        if protein_peptide_index is not None
        else peptide_bed_dict
    )
    for peptide_name in candidate_peptides:
        locations = peptide_bed_dict.get(peptide_name, [])
        if not locations or not any(
            _location_matches_protein_and_model(location, protein_id, pred_info)
            for location in locations
        ):
            continue

        def supports_same_terminal(location):
            blocks = location.get('blocks', [])
            if (
                not blocks
                or location.get('seqid') != pred_info.get('seqid')
                or location.get('strand') != strand
            ):
                return False

            if terminal == 'start':
                # pep_start=1 means there is no preceding amino acid and thus no
                # preceding R/K tryptic cleavage site.
                if (
                    location.get('pep_start') != 1
                    or location.get('previous_aa', '').upper() in {'R', 'K'}
                ):
                    return False
                peptide_boundary = min(start for start, _ in blocks) if strand == '+' else max(end for _, end in blocks)
                return peptide_boundary == target

            if location.get('pep_end') != location.get('prot_len'):
                return False
            if peptide_name.upper().endswith(('R', 'K')):
                return False
            peptide_boundary = max(end for _, end in blocks) if strand == '+' else min(start for start, _ in blocks)
            # A terminal stop codon is inconsistently included as three CDS
            # bases, so a terminal peptide may end at the boundary or 3 nt in.
            return abs(peptide_boundary - target) <= 3

        # A peptide mapping internally or to another genomic event has an
        # ambiguous origin and cannot support this event. The same peptide may
        # map to multiple isoforms only when every mapping supports this exact
        # genomic terminus.
        if all(supports_same_terminal(location) for location in locations):
            supporting.add(peptide_name)

    return supporting


def get_peptides_supporting_novel_splice_sites(
    protein_id, pred_info, reference_event_index, peptide_bed_dict,
    protein_peptide_index=None,
):
    """Return peptide-supported junctions with novel donor and acceptor sites."""
    if not pred_info:
        return {}
    key = (pred_info.get('seqid'), pred_info.get('strand'))
    reference_sites = reference_event_index.get(key, {})
    reference_donors = reference_sites.get('donors', set())
    reference_acceptors = reference_sites.get('acceptors', set())
    strand = pred_info.get('strand')
    supported = defaultdict(set)
    candidate_peptides = (
        protein_peptide_index.get(protein_id, set())
        if protein_peptide_index is not None
        else peptide_bed_dict
    )

    for junction in get_splice_junctions(pred_info.get('cds_regions', [])):
        donor, acceptor = _donor_acceptor(junction, strand)
        if donor in reference_donors or acceptor in reference_acceptors:
            continue
        for peptide_name in candidate_peptides:
            locations = peptide_bed_dict.get(peptide_name, [])
            supports_current_protein = any(
                _location_matches_protein_and_model(location, protein_id, pred_info)
                and peptide_spans_junction(location, junction)
                for location in locations
            )
            supports_only_this_event = bool(locations) and all(
                location.get('seqid') == pred_info.get('seqid')
                and location.get('strand') == strand
                and peptide_spans_junction(location, junction)
                for location in locations
            )
            if supports_current_protein and supports_only_this_event:
                supported[junction].add(peptide_name)
    return dict(supported)


def load_peptides_bed(bed_file):
    """
    Load peptides_mapped.bed file and return peptide locations.
    
    Returns: dict mapping peptide_name to list of genomic locations
        {
            peptide_name: [{
                'seqid': chr,
                'start': pos,
                'end': pos,
                'strand': strand,
                'protein': protein_id,
                'blocks': [(block_start_1based, block_end_1based), ...]
            }, ...]
        }
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
            name = fields[3]  # peptide sequence in N-to-C orientation
            score = fields[4]
            strand = fields[5]

            # Parse BED12 blocks when available.
            # BED is 0-based, end-exclusive. Convert blocks to 1-based inclusive
            # to match GTF-derived junction coordinates.
            blocks = []
            if len(fields) >= 12:
                try:
                    block_count = int(fields[9])
                    block_sizes = [int(x) for x in fields[10].rstrip(',').split(',') if x]
                    block_starts = [int(x) for x in fields[11].rstrip(',').split(',') if x]
                    if block_count == len(block_sizes) == len(block_starts):
                        for rel_start, size in zip(block_starts, block_sizes):
                            b_start_0 = start + rel_start
                            b_end_0 = b_start_0 + size
                            blocks.append((b_start_0 + 1, b_end_0))
                except ValueError:
                    blocks = []
            
            metadata = {}
            for attr in fields[12:]:
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    metadata[key] = value

            def _metadata_int(key):
                try:
                    return int(metadata[key])
                except (KeyError, TypeError, ValueError):
                    return None
            
            peptide_locations[name].append({
                'seqid': seqid,
                'start': start,
                'end': end,
                'strand': strand,
                'protein': metadata.get('protein'),
                'blocks': blocks,
                'pep_start': _metadata_int('pep_start'),
                'pep_end': _metadata_int('pep_end'),
                'prot_len': _metadata_int('prot_len'),
                'previous_aa': metadata.get('previous_aa', ''),
            })
    
    return peptide_locations


def load_genomic_model_evidence(genomic_bed, projection_bed):
    """Attach model metadata only to matching genomic placements.

    Retain every additional genomic placement when testing global ambiguity.
    Historical protein normalization uses I->L, so join on that equivalence class.
    Unsupported protein projections must not become genomic event evidence.
    """
    genomic = load_peptides_bed(genomic_bed)
    projected = load_peptides_bed(projection_bed)
    metadata = defaultdict(list)
    def key(name, location):
        return (name.replace('I', 'L'), location['seqid'], location['strand'],
                tuple(map(tuple, location['blocks'])))
    for name, locations in projected.items():
        for location in locations:
            metadata[key(name, location)].append(location)
    merged = defaultdict(list)
    seen = set()
    for name, locations in genomic.items():
        for location in locations:
            token = key(name, location)
            if token in seen:
                continue
            seen.add(token)
            matches = metadata.get(token)
            merged[token[0]].extend(matches if matches else [location])
    return merged


def peptide_spans_junction(peptide_loc, junction):
    """
    Check if a peptide spans a specific splice junction.
    
    Args:
        peptide_loc: dict with 'blocks' parsed from BED12
        junction: tuple (exon_end, next_exon_start)

    Returns: True only when peptide has a split mapping with adjacent blocks that
    exactly flank the junction boundaries.
    """
    exon_end, next_exon_start = junction
    blocks = peptide_loc.get('blocks', [])

    # Without split blocks, a peptide cannot provide true splice-junction support.
    if len(blocks) < 2:
        return False

    blocks = sorted(blocks, key=lambda x: x[0])
    for i in range(len(blocks) - 1):
        left_end = blocks[i][1]
        right_start = blocks[i + 1][0]
        if left_end == exon_end and right_start == next_exon_start:
            return True

    return False


def get_corresponding_ref_junctions(pred_junction, ref_only_junctions):
    """
    Get reference-only junctions that correspond to a predicted-only junction.

    "Corresponding" means they share one splice boundary (same donor or acceptor)
    and differ on the other boundary.
    """
    pred_exon_end, pred_next_exon_start = pred_junction
    corresponding = []

    for ref_junction in ref_only_junctions:
        ref_exon_end, ref_next_exon_start = ref_junction
        if pred_exon_end == ref_exon_end or pred_next_exon_start == ref_next_exon_start:
            corresponding.append(ref_junction)

    return corresponding


def get_peptides_supporting_different_junctions(protein_id, pred_info, ref_info, peptide_bed_dict):
    """
    Find peptides that discriminatively support predicted-only splice junctions.

    A peptide is counted only if it:
    1) spans a predicted-only junction, and
    2) does NOT span the corresponding reference-only junction(s).
    
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
    
    # Predicted-only and reference-only junction sets
    different_junctions = pred_junctions - ref_junctions
    ref_only_junctions = ref_junctions - pred_junctions
    
    if not different_junctions:
        return 0
    
    # Some BED files do not carry protein IDs. In that case, fall back to
    # chromosome-only matching instead of returning zero for every query.
    has_any_protein_ids = any(
        loc.get('protein')
        for locations in peptide_bed_dict.values()
        for loc in locations
    )

    # Map each predicted-only junction to corresponding reference-only junction(s)
    corresponding_ref_by_pred = {
        junction: get_corresponding_ref_junctions(junction, ref_only_junctions)
        for junction in different_junctions
    }

    # Find discriminative peptides for this protein
    supporting_peptides = set()
    
    # Search through all peptides in the BED file
    for peptide_name, locations in peptide_bed_dict.items():
        for loc in locations:
            # Only count peptides mapped to this protein on the same chromosome.
            if loc['seqid'] != seqid or loc.get('strand') != pred_info.get('strand'):
                continue

            if has_any_protein_ids and loc.get('protein') not in (None, protein_id):
                continue

            # Count only discriminative support: spans pred-only junction and not
            # the corresponding ref-only junction(s).
            for junction in different_junctions:
                if not peptide_spans_junction(loc, junction):
                    continue

                corresponding_ref = corresponding_ref_by_pred.get(junction, [])
                if any(peptide_spans_junction(loc, ref_j) for ref_j in corresponding_ref):
                    continue

                supporting_peptides.add(peptide_name)
                break
    
    return len(supporting_peptides)


def generate_annotation_report(output_dir, supported_gtf, reference_gtf, 
                               tmap_file=None, all_proteins_scores=None, peptides_bed=None,
                               protein_fasta=None):
    """
    Generate comprehensive annotation comparison report.
    
    Report event-specific peptide support at coding boundaries absent from every
    reference isoform. Novel-gene export is handled separately by strict CDS
    non-overlap in parse_gffcompare_tmap.py.
    
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
    compare_dir = output_dir / "Compare_to_Reference"

    # Auto-detect file paths if not provided
    if tmap_file is None:
        tmap_candidates = (
            [
                compare_dir / "gffcmp.peptide_mapped_proteins_cleaned.gtf.tmap",
                output_dir / "gffcmp.peptide_mapped_proteins_cleaned.gtf.tmap",
                compare_dir / "gffcmp.supported_proteins_cleaned.gtf.tmap",
                output_dir / "gffcmp.supported_proteins_cleaned.gtf.tmap",
                compare_dir / "gffcmp.supported_proteins.gtf.tmap",
                output_dir / "gffcmp.supported_proteins.gtf.tmap",
                output_dir / "gffcompare.supported_proteins.gtf.tmap",
            ]
            + sorted(compare_dir.glob("*.tmap"))
            + sorted(output_dir.glob("*.tmap"))
        )
    else:
        provided_tmap = Path(tmap_file)
        tmap_candidates = [
            provided_tmap,
            compare_dir / provided_tmap.name,
            output_dir / provided_tmap.name,
        ]

    tmap_file = _pick_tmap_file(tmap_candidates) or Path(tmap_candidates[0])

    if all_proteins_scores is None:
        all_proteins_scores = output_dir / "all_proteins_scores.tsv"

    if peptides_bed is None:
        peptides_bed = output_dir / "peptides_mapped.bed"

    all_proteins_scores = Path(all_proteins_scores)
    peptides_bed = Path(peptides_bed)
    supported_gtf = Path(supported_gtf)
    reference_gtf = Path(reference_gtf)

    if not tmap_file.exists():
        looked_in = "\n  - ".join(str(path) for path in tmap_candidates)
        raise FileNotFoundError(f"TMAP file not found. Looked in:\n  - {looked_in}")

    print(f"Loading TMAP file: {tmap_file}")
    tmap_df = pd.read_csv(tmap_file, sep='\t')
    
    print(f"Parsing GTF files...")
    supported_cds = parse_gtf_to_cds_dict(supported_gtf)
    reference_cds = parse_gtf_to_cds_dict(reference_gtf)
    reference_event_index = build_reference_event_index(reference_cds)
    
    print(f"Loading peptides BED file: {peptides_bed}")
    projection_bed = output_dir / "protein_projection" / "peptides_mapped.bed"
    if projection_bed.exists() and peptides_bed == output_dir / "peptides_mapped.bed":
        peptide_bed_dict = load_genomic_model_evidence(peptides_bed, projection_bed)
    else:
        peptide_bed_dict = load_peptides_bed(peptides_bed)
    protein_peptide_index = build_protein_peptide_index(peptide_bed_dict)
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
    
    # Event discovery is peptide-driven and must not be gated by the protein
    # classifier. Iterate over every peptide-mapped model in the score table;
    # gffcompare metadata is attached when available but does not determine
    # eligibility for an event call.
    tmap_by_query = {
        str(row.qry_id): row
        for row in tmap_df.drop_duplicates(subset=['qry_id']).itertuples()
    }
    candidate_proteins = {
        str(protein_id)
        for protein_id, score_values in scores_dict.items()
        if pd.to_numeric(score_values.get('mapped_peptides', 0), errors='coerce') > 0
    }
    if not candidate_proteins:
        candidate_proteins = set(protein_peptide_index)

    for qry_id in sorted(candidate_proteins):
        tmap_row = tmap_by_query.get(qry_id)
        ref_id = getattr(tmap_row, 'ref_id', '-')
        class_code = getattr(tmap_row, 'class_code', 'not_evaluated')
        
        pred_info = supported_cds.get(qry_id)
        protein_id = pred_info.get('protein_id') if pred_info else qry_id
        if protein_id not in scores_dict and qry_id in scores_dict:
            protein_id = qry_id

        # Get peptide support metrics
        if protein_id in scores_dict:
            peptides = scores_dict[protein_id].get('mapped_peptides', 0)
            coverage = scores_dict[protein_id].get('sequence_coverage', 0)
            n_term_peptides = scores_dict[protein_id].get('N_terminal_peptides', 0)
            c_term_peptides = scores_dict[protein_id].get('C_terminal_peptides', 0)
            splice_peptides = scores_dict[protein_id].get('splice_peptides', 0)
        else:
            peptides = 0
            coverage = 0
            n_term_peptides = 0
            c_term_peptides = 0
            splice_peptides = 0
        
        # Get predicted CDS info. A prediction without
        # same-strand CDS overlap is handled as a novel/antisense candidate,
        # not as an alternative event at a known coding locus.
        if not pred_info or not _has_same_strand_reference_cds_overlap(pred_info, reference_event_index):
            continue

        key = (pred_info.get('seqid'), pred_info.get('strand'))
        reference_boundaries = reference_event_index.get(key, {})
        pred_start, pred_stop = get_first_and_last_cds(
            pred_info.get('cds_regions', []), pred_info.get('strand')
        )

        events = []

        start_support = get_terminal_supporting_peptides(
            protein_id, pred_info, peptide_bed_dict, 'start',
            protein_peptide_index,
        )
        if pred_start not in reference_boundaries.get('starts', set()) and start_support:
            events.append((
                'peptide_support_different_start',
                f"start:{pred_start}",
                start_support,
            ))

        stop_support = get_terminal_supporting_peptides(
            protein_id, pred_info, peptide_bed_dict, 'stop',
            protein_peptide_index,
        )
        reference_stops = reference_boundaries.get('stops', set())
        stop_absent_from_reference = not any(
            abs(pred_stop - reference_stop) <= 3
            for reference_stop in reference_stops
        )
        if stop_absent_from_reference and stop_support:
            events.append((
                'peptide_support_different_stop',
                f"stop:{pred_stop}",
                stop_support,
            ))

        supported_junctions = get_peptides_supporting_novel_splice_sites(
            protein_id, pred_info, reference_event_index, peptide_bed_dict,
            protein_peptide_index,
        )
        if supported_junctions:
            junction_peptides = set().union(*supported_junctions.values())
            junction_text = ';'.join(
                f"junction:{left_end}-{right_start}"
                for left_end, right_start in sorted(supported_junctions)
            )
            events.append((
                'peptide_support_different_splice',
                junction_text,
                junction_peptides,
            ))

        for category, event_coordinates, evidence_peptides in events:
            category_counts[category]['count'] += 1
            category_counts[category]['ids'].add(protein_id)
            report_data.append({
                'predicted_id': protein_id,
                'gffcompare_query_id': qry_id,
                'reference_id': ref_id,
                'gffcompare_class_code': class_code,
                'structural_difference': category.removeprefix('peptide_support_'),
                'peptide_evidence_category': category,
                'event_coordinates': event_coordinates,
                'event_supporting_peptides': ','.join(sorted(evidence_peptides)),
                'N_terminal_peptides': len(start_support),
                'C_terminal_peptides': len(stop_support),
                'splice_peptides': len(set().union(*supported_junctions.values())) if supported_junctions else 0,
                'mapped_peptides': peptides,
                'sequence_coverage': coverage
            })
    
    # Generate TSV report
    report_columns = [
        'predicted_id', 'gffcompare_query_id', 'reference_id',
        'gffcompare_class_code', 'structural_difference',
        'peptide_evidence_category', 'event_coordinates',
        'event_supporting_peptides', 'N_terminal_peptides',
        'C_terminal_peptides', 'splice_peptides', 'mapped_peptides',
        'sequence_coverage',
    ]
    report_df = pd.DataFrame(report_data, columns=report_columns)
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
    
    summary_df = pd.DataFrame(
        summary_report,
        columns=['Category', 'Count', 'Protein_IDs'],
    )
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"Summary report written to: {summary_file}")
    
    # Extract GTF entries by category
    print(f"\nExtracting GTF entries by category...")
    gtf_dir = compare_dir
    different_dir = compare_dir / "Different"
    different_dir.mkdir(parents=True, exist_ok=True)

    # Remove stale category files from previous runs.
    for stale in different_dir.glob("*.gtf"):
        stale.unlink()
    for stale in different_dir.glob("*.faa"):
        stale.unlink()
    
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
    
    # Load protein sequences from the provided protein FASTA (or auto-detect).
    protein_sequences = {}
    if not protein_fasta:
        autodetected_fasta = output_dir / "supported_proteins.fa"
        if autodetected_fasta.exists():
            protein_fasta = str(autodetected_fasta)
        else:
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
            # Extract protein/transcript ID from either GFF3-style or GTF-style attributes.
            attrs = entry['attributes']
            candidates = (
                extract_attribute(attrs, 'protein_id'),
                extract_attribute(attrs, 'transcript_id'),
                extract_attribute(attrs, 'Parent'),
                extract_attribute(attrs, 'ID'),
            )
            protein_id = next(
                (candidate for candidate in candidates if candidate in protein_ids),
                None,
            )

            if protein_id in protein_ids:
                category_gtf_lines.append(entry['line'])
        
        # Write GTF file for this category
        if category_gtf_lines:
            # Sanitize filename
            safe_category = category.replace('/', '_').replace(' ', '_')
            category_dir = different_dir
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
