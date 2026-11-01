"""Anchored, exactly-one-amino-acid-edit rescue of unmapped peptides.

A query is partitioned into two non-overlapping exact seeds. At edit distance one,
at least one seed survives, including when the edit is an indel. Candidate starts
are expanded by +/-1 amino acid and verified against target lengths L-1, L, L+1.
All standard-amino-acid matches are retained; stops/unresolved codons are barriers.
"""
from collections import Counter, defaultdict

import ahocorasick
import pysam

from .io import AMINO_ACIDS
from .sequence import canonical, project, reverse_complement, translate

_WORKER = None


def one_edit(query, target, il_equivalent=True):
    """Return deterministic full-query alignment details for distance exactly 1.

    CIGAR: '=' consumes both sequences, X substitutes, I consumes query only,
    D consumes assembly translation only. Positions are 1-based in the query;
    a terminal D is positioned just after the query. I/L equivalence costs zero.
    """
    q, t = canonical(query, il_equivalent), canonical(target, il_equivalent)
    delta = len(q) - len(t)
    if abs(delta) > 1 or not target or set(target) - AMINO_ACIDS:
        return None
    position = 0
    bound = min(len(q), len(t))
    while position < bound and q[position] == t[position]:
        position += 1
    if delta == 0:
        if position == len(q) or q[position+1:] != t[position+1:]:
            return None
        operation, kind = "X", "substitution"
        query_aa, reference_aa = query[position], target[position]
        tail = len(q) - position - 1
    elif delta == 1:
        if q[position+1:] != t[position:]:
            return None
        operation, kind = "I", "insertion"
        query_aa, reference_aa = query[position], "-"
        tail = len(q) - position - 1
    else:
        if q[position:] != t[position+1:]:
            return None
        operation, kind = "D", "deletion"
        query_aa, reference_aa = "-", target[position]
        tail = len(q) - position
    cigar = (f"{position}=" if position else "") + "1" + operation + (f"{tail}=" if tail else "")
    return {"edit_distance": 1, "edit_type": kind, "edit_query_position": position+1,
            "query_aa": query_aa, "reference_aa": reference_aa, "alignment_cigar": cigar}


def make_rescue_automaton(peptides, peptide_ids, il_equivalent=True):
    grouped = defaultdict(set)
    for pid in peptide_ids:
        query = canonical(peptides[pid], il_equivalent)
        midpoint = len(query) // 2
        if midpoint < 1:
            raise ValueError("one-edit rescue requires peptides of at least two amino acids")
        grouped[query[:midpoint]].add((pid, 0))
        grouped[query[midpoint:]].add((pid, midpoint))
    if not grouped:
        return None
    automaton = ahocorasick.Automaton()
    for seed, positions in grouped.items():
        automaton.add_word(seed, (len(seed), tuple(sorted(positions))))
    automaton.make_automaton()
    return automaton


def protein_candidates(protein, automaton, peptides, il_equivalent=True, stats=None):
    """Yield (query_id, target_start_aa, matched_sequence, edit_details)."""
    if automaton is None:
        return
    stats = stats if stats is not None else Counter()
    searchable = canonical(protein, il_equivalent)
    seen = {}
    stride = len(protein) + 1
    for end, (seed_length, positions) in automaton.iter(searchable):
        stats["seed_occurrences"] += 1
        seed_start = end - seed_length + 1
        for pid, offset in positions:
            query = peptides[pid]
            size = len(query)
            base = seed_start - offset
            # Only propose alignments where this seed is the untouched half.
            # Prefix seed: same start for all edit types. Suffix seed: the
            # start shifts with the indel so the seed still ends at query end.
            # Every distance-1 alignment has at least one such untouched half.
            proposals = ((base, size, 1),
                         (base if offset == 0 else base+1, size-1, 2),
                         (base if offset == 0 else base-1, size+1, 4))
            for start, target_length, bit in proposals:
                if start < 0 or start + target_length > len(protein):
                    continue
                token = pid * stride + start
                mask = seen.get(token, 0)
                if mask & bit:
                    continue
                seen[token] = mask | bit
                if not mask:
                    stats["candidate_starts"] += 1
                target = protein[start:start+target_length]
                stats["candidate_alignments"] += 1
                details = one_edit(query, target, il_equivalent)
                if details is not None:
                    stats["verified_alignments_before_deduplication"] += 1
                    yield pid, start, target, details


def search_rescue_dna(dna, chrom, strand, segments, automaton, peptides, table=1,
                      il_equivalent=True, required_junctions=0, support=0,
                      motifs="", source="genome", stats=None):
    boundaries = []
    cumulative = 0
    for a,b in segments[:-1]:
        cumulative += b-a
        boundaries.append(cumulative)
    for frame in range(3):
        protein = translate(dna[frame:], table)
        for pid, start, matched, details in protein_candidates(protein, automaton, peptides, il_equivalent, stats):
            nt_start, nt_end = frame+3*start, frame+3*(start+len(matched))
            if required_junctions and not all(nt_start < boundary < nt_end for boundary in boundaries):
                continue
            blocks = project(nt_start, nt_end, segments, strand)
            yield pid, chrom, strand, blocks, matched, source, support, motifs, details


def init_rescue_worker(genome, index, automaton, peptides, table, il_equivalent):
    global _WORKER
    _WORKER = (pysam.FastaFile(genome, filepath_index=index), automaton, peptides, table, il_equivalent)


def rescue_chunk(task):
    fasta, automaton, peptides, table, il_equivalent = _WORKER
    chrom, start, core_end, fetch_end = task
    forward = fasta.fetch(chrom, start, fetch_end).upper()
    stats, hits = Counter(), []
    for strand, dna in (("+", forward), ("-", reverse_complement(forward))):
        for hit in search_rescue_dna(dna, chrom, strand, ((start,fetch_end),), automaton, peptides, table, il_equivalent, stats=stats):
            if start <= hit[3][0][0] < core_end:
                hits.append(hit)
    return hits, dict(stats)


def rescue_path(fasta, path, automaton, peptides, table=1, il_equivalent=True, stats=None):
    chrom, strand, segments, support, motifs, source = path
    if any(a>=b for a,b in segments):
        return
    dna = "".join(fasta.fetch(chrom,a,b) if strand=="+" else reverse_complement(fasta.fetch(chrom,a,b)) for a,b in segments).upper()
    yield from search_rescue_dna(dna, chrom, strand, segments, automaton, peptides, table,
                                il_equivalent, len(segments)-1, support, motifs, source, stats)
