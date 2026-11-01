"""Aho-Corasick searches of chunked six-frame translations."""
from collections import defaultdict

import ahocorasick
import pysam

from .sequence import canonical, project, reverse_complement, translate

_WORKER = None


def make_automaton(peptides, il_equivalent=True):
    grouped = defaultdict(list)
    for index, peptide in enumerate(peptides):
        grouped[canonical(peptide, il_equivalent)].append(index)
    automaton = ahocorasick.Automaton()
    for key, indexes in grouped.items():
        automaton.add_word(key, (len(key), tuple(indexes)))
    automaton.make_automaton()
    return automaton


def init_worker(genome, index, automaton, peptides, table, il_equivalent):
    global _WORKER
    _WORKER = (pysam.FastaFile(genome, filepath_index=index), automaton, peptides, table, il_equivalent)


def search_dna(dna, chrom, strand, segments, automaton, peptides, table=1, il_equivalent=True, required_junctions=0, support=0, motifs="", source="genome"):
    boundaries = []
    total = 0
    for a, b in segments[:-1]:
        total += b - a
        boundaries.append(total)
    for frame in range(3):
        protein = translate(dna[frame:], table)
        searchable = canonical(protein, il_equivalent)
        for end, (size, indexes) in automaton.iter(searchable):
            aa_start = end - size + 1
            nt_start, nt_end = frame + 3 * aa_start, frame + 3 * (end + 1)
            # A path must cross every supplied edge; subpaths are searched separately.
            if required_junctions and not all(nt_start < b < nt_end for b in boundaries):
                continue
            blocks = project(nt_start, nt_end, segments, strand)
            matched = protein[aa_start:end + 1]
            for pid in indexes:
                yield (pid, chrom, strand, blocks, matched, source, support, motifs)


def search_chunk(task):
    fasta, automaton, peptides, table, il_equivalent = _WORKER
    chrom, start, core_end, fetch_end = task
    forward = fasta.fetch(chrom, start, fetch_end).upper()
    hits = []
    for strand, dna in (("+", forward), ("-", reverse_complement(forward))):
        for hit in search_dna(dna, chrom, strand, ((start, fetch_end),), automaton, peptides, table, il_equivalent):
            # Own each genomic lower coordinate exactly once, on either strand.
            if start <= hit[3][0][0] < core_end:
                hits.append(hit)
    return hits


def chunk_tasks(lengths, chunk_size, max_length):
    overlap = 3 * max_length - 1
    for chrom, size in lengths.items():
        for start in range(0, size, chunk_size):
            end = min(size, start + chunk_size)
            yield chrom, start, end, min(size, end + overlap)


def search_path(fasta, path, automaton, peptides, table=1, il_equivalent=True):
    chrom, strand, segments, support, motifs, source = path
    if any(a >= b for a, b in segments):
        return
    dna = "".join(fasta.fetch(chrom, a, b) if strand == "+" else reverse_complement(fasta.fetch(chrom, a, b)) for a, b in segments).upper()
    yield from search_dna(dna, chrom, strand, segments, automaton, peptides, table, il_equivalent, len(segments) - 1, support, motifs, source)
