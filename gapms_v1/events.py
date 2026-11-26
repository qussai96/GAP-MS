"""Shared evidence rules for annotation-change events.

Every event type (alternative start/stop/splice, terminal extensions, novel
ORFs) is judged by the same peptide rules:

* a peptide is *diagnostic* only if no supplied annotation explains it at the
  same genomic locus;
* a diagnostic peptide counts as specific only if it has exactly one genomic
  locus among **all** exact placements (tryptic or not);
* specific peptides with failed identification q-values do not support a
  change. Peptide length is reported, not filtered: whether a short unique
  peptide is a chance match is left to the user, aided by the
  length-stratified decoy FDR of the local search;
* event-specific blockers (e.g. conflicting upstream evidence) override support.

Rows describing the same biological feature are merged into one row.
"""
from collections import defaultdict
import json

from peplocus.sequence import canonical
from .common import stable_id

SUPPORTED = ('high_confidence_candidate', 'rna_peptide_supported', 'peptide_supported')
RANK = {'high_confidence_candidate': 4, 'rna_peptide_supported': 3, 'peptide_supported': 2, 'ambiguous': 1, 'not_supported': 0}
# A blocker found for any contributing model applies to the merged feature.
BLOCKERS = {
    'conflicting_upstream_evidence', 'initiator_methionine_ambiguity', 'upstream_stop_without_start',
    'noncanonical_intron_in_candidate', 'candidate_sequence_or_phase_issues', 'trypsin_site_before_proposed_start',
}
# Not a hard blocker: another contributing row may hold a path-specific peptide.
def quality_lookup(metadata):
    qualities = defaultdict(set)
    for peptide, item in metadata.items():
        qualities[canonical(peptide)].add(item['quality'])
    return qualities


def _junctions(blocks):
    blocks=sorted(tuple(b) for b in blocks)
    return {(a[1], b[0]) for a, b in zip(blocks, blocks[1:]) if a[1] < b[0]}


def _five_prime(blocks, strand):
    blocks=sorted(tuple(b) for b in blocks)
    return blocks[0][0] if strand == '+' else blocks[-1][1]


def _three_prime(blocks, strand):
    blocks=sorted(tuple(b) for b in blocks)
    return blocks[-1][1] if strand == '+' else blocks[0][0]


def feature_consistent(kind, subtype, strand, placement, members, junctions=(), host=None):
    """Does every alternative path of the peptide at this locus imply the feature?

    A peptide that two different same-locus paths explain (e.g. two donors
    that translate identically) proves that one of them exists, not which.
    """
    if kind == 'alternative_splice' and subtype in {'new_junction', 'junction_shift', 'exon_skip'}:
        wanted = {tuple(j) for j in junctions}
        return all(wanted <= _junctions(m) for m in members)
    if kind == 'alternative_splice' and subtype in {'retained_intron', 'novel_exon'} and host:
        return all(any(host[0] <= a and b <= host[1] for a, b in m) for m in members)
    if kind in {'alternative_start', 'n_terminal_extension'}:
        return all(_five_prime(m, strand) == _five_prime(placement, strand) for m in members)
    if kind in {'alternative_stop', 'c_terminal_extension'}:
        return all(_three_prime(m, strand) == _three_prime(placement, strand) for m in members)
    return True


def judge(sequences, strict_counts, qualities, blockers=(), consistent=None):
    """Return (base_confidence, specific_sequences, reasons).

    base_confidence is 'peptide_supported', 'ambiguous' or 'not_supported';
    callers may promote supported rows using RNA or structural evidence.
    """
    sequences = set(sequences)
    unique = {p for p in sequences if strict_counts.get(p, 0) == 1}
    reasons = []
    specific = set()
    for p in unique:
        quality = qualities.get(p, {'unknown'})
        if quality == {'qvalue_fail'}:
            continue
        specific.add(p)
    path_ambiguous = set()
    if consistent is not None:
        path_ambiguous = {p for p in specific if p not in consistent}
        specific -= path_ambiguous
    if not unique:
        reasons.append('no_single_locus_diagnostic_peptide')
        confidence = 'ambiguous'
    elif not specific:
        if any(qualities.get(p) == {'qvalue_fail'} for p in unique):
            reasons.append('single_locus_peptides_failed_identification_qvalue')
        if path_ambiguous:
            reasons.append('peptide_path_ambiguous_within_locus')
        confidence = 'not_supported'
    else:
        confidence = 'peptide_supported'
    blockers = [b for b in blockers if b]
    if blockers:
        reasons.extend(sorted(set(blockers)))
        confidence = 'not_supported'
    return confidence, specific, reasons


def feature_id(kind, chrom, strand, feature):
    return stable_id('F', kind, chrom, strand, feature)


def _as_list(value):
    if value in (None, ''):
        return []
    if isinstance(value, str):
        return json.loads(value)
    return list(value)


def merge_events(rows):
    """Collapse rows with the same event_id (one biological feature).

    The best-supported contributing row provides the confidence, unless any
    contributing row carries a hard blocker. Evidence lists are unioned.
    """
    grouped = defaultdict(list)
    for row in rows:
        grouped[row['event_id']].append(row)
    merged = []
    for event_id, items in grouped.items():
        items = sorted(items, key=lambda r: (-RANK.get(r['confidence'], 0), r.get('evidence_origin', ''), str(r.get('candidate_model_id', ''))))
        base = dict(items[0])
        if len(items) > 1:
            def union(key):
                values = []
                for r in items:
                    for v in _as_list(r.get(key)):
                        if v not in values:
                            values.append(v)
                return values
            for key in ['peptides', 'mapping_ids', 'reference_model_ids', 'annotation_sources', 'limitations',
                        'whole_model_review_reasons', 'confidence_reasons', 'supporting_model_ids']:
                base[key] = sorted(union(key), key=str)
            base['local_peptide_blocks'] = sorted({tuple(b) for r in items for b in _as_list(r.get('local_peptide_blocks'))})
            base['diagnostic_junctions'] = sorted({tuple(b) for r in items for b in _as_list(r.get('diagnostic_junctions'))})
            junctions = {}
            for r in items:
                for j in _as_list(r.get('rna_junctions')):
                    junctions[(j['chrom'], j['start'], j['end'], j['strand'])] = j
            base['rna_junctions'] = [junctions[k] for k in sorted(junctions)]
            for key in ['mapping_multiplicity', 'peptide_loci', 'reference_cds_blocks']:
                combined = {}
                for r in items:
                    value = r.get(key)
                    combined.update(json.loads(value) if isinstance(value, str) and value else (value or {}))
                base[key] = combined
            contexts = []
            for r in items:
                for c in _as_list(r.get('cleavage_context')):
                    if c not in contexts:
                        contexts.append(c)
            base['cleavage_context'] = contexts
            base['evidence_origin'] = ','.join(sorted({o for r in items for o in str(r.get('evidence_origin', '')).split(',') if o}))
            base['distinct_diagnostic_peptides'] = len({canonical(p) for p in base['peptides']})
            upstream = [r.get('upstream_peptide_support') for r in items if r.get('upstream_peptide_support') not in (None, '')]
            base['upstream_peptide_support'] = max(int(u) for u in upstream) if upstream else None
            blocked = sorted(BLOCKERS & set(base['confidence_reasons']))
            if blocked:
                base['confidence'] = 'not_supported'
            structures = {json.dumps(_as_list(r.get('candidate_cds_blocks'))) for r in items if _as_list(r.get('candidate_cds_blocks'))}
            if len(structures) > 1:
                base['limitations'] = sorted(set(base['limitations']) | {'multiple_candidate_coding_structures_share_this_evidence'})
                if base['confidence'] == 'high_confidence_candidate':
                    base['confidence'] = 'peptide_supported'
                    base['confidence_reasons'] = sorted(set(base['confidence_reasons']) | {'coding_structure_not_resolved_between_isoforms'})
        merged.append(base)
    return merged
