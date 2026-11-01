"""RNA-seq junction extraction and bounded splice-path construction."""
from bisect import bisect_left
from collections import Counter, defaultdict
from dataclasses import dataclass
import csv

import pysam

from .sequence import reverse_complement

CANONICAL = {"GT-AG", "GC-AG", "AT-AC"}


@dataclass(frozen=True, order=True)
class Junction:
    chrom: str
    start: int
    end: int
    strand: str
    support: int = 0
    source: str = "bam"
    motif: str = ""


def motif(fasta, chrom, start, end, strand):
    if strand == "+":
        return fasta.fetch(chrom, start, start + 2).upper() + "-" + fasta.fetch(chrom, end - 2, end).upper()
    return reverse_complement(fasta.fetch(chrom, end - 2, end).upper()) + "-" + reverse_complement(fasta.fetch(chrom, start, start + 2).upper())


def resolve_junction(fasta, chrom, start, end, strand, support, source):
    if strand in {"+", "-"}:
        return [Junction(chrom, start, end, strand, support, source, motif(fasta, chrom, start, end, strand))]
    candidates = [Junction(chrom, start, end, s, support, source, motif(fasta, chrom, start, end, s)) for s in ("+", "-")]
    canonical = [j for j in candidates if j.motif in CANONICAL]
    return canonical or candidates


def extract_bam(path, fasta, lengths, min_mapq=10, min_anchor=8, min_support=2, min_intron=20, max_intron=500000):
    counts = Counter()
    stats = Counter()
    with pysam.AlignmentFile(str(path), "rb") as bam:
        common = set(bam.references) & lengths.keys()
        if not common:
            raise ValueError("BAM and assembly have no matching sequence names")
        for chrom, length in zip(bam.references, bam.lengths):
            if chrom in lengths and lengths[chrom] != length:
                raise ValueError(f"BAM/assembly length mismatch for {chrom}: {length} versus {lengths[chrom]}")
        for read in bam.fetch(until_eof=True):
            stats["alignments_seen"] += 1
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail or read.is_duplicate or read.mapping_quality < min_mapq:
                continue
            cigar = read.cigartuples
            if not cigar or not any(op == 3 for op, _ in cigar):
                continue
            chrom = read.reference_name
            if chrom not in lengths:
                stats["spliced_alignments_unknown_contig"] += 1
                continue
            strand = "."
            for tag in ("XS", "ts"):
                if read.has_tag(tag) and read.get_tag(tag) in {"+", "-"}:
                    strand = read.get_tag(tag)
                    break
            position = read.reference_start
            for i, (op, size) in enumerate(cigar):
                if op == 3:
                    left = right = 0
                    k = i - 1
                    while k >= 0 and cigar[k][0] in {0, 7, 8}:
                        left += cigar[k][1]
                        k -= 1
                    k = i + 1
                    while k < len(cigar) and cigar[k][0] in {0, 7, 8}:
                        right += cigar[k][1]
                        k += 1
                    if min_intron <= size <= max_intron and min(left, right) >= min_anchor:
                        counts[(chrom, position, position + size, strand)] += 1
                        stats["junction_observations"] += 1
                if op in {0, 2, 3, 7, 8}:
                    position += size
    merged = Counter()
    # Resolve before applying minimum support so tagged and untagged reads combine.
    for (chrom, start, end, strand), count in counts.items():
        for j in resolve_junction(fasta, chrom, start, end, strand, count, "bam"):
            merged[(j.chrom, j.start, j.end, j.strand, j.motif)] += count
    junctions = [Junction(chrom, start, end, strand, count, "bam", m) for (chrom, start, end, strand, m), count in merged.items() if count >= min_support]
    stats["junctions_retained"] = len(junctions)
    return sorted(junctions), dict(stats)


def read_junctions(path, fasta, lengths, min_support=1):
    result = []
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chrom", "start", "end", "strand", "support"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError("junction TSV requires chrom, start, end, strand, support columns (0-based intron boundaries)")
        for row in reader:
            chrom, start, end = row["chrom"], int(row["start"]), int(row["end"])
            if chrom not in lengths or not (0 <= start < end <= lengths[chrom]):
                raise ValueError(f"junction outside assembly: {row}")
            if row["strand"] not in {"+", "-", "."}:
                raise ValueError(f"invalid junction strand: {row['strand']}")
            count = int(row["support"])
            if count >= min_support or "annotation" in row.get("source", "").split(","):
                result.extend(resolve_junction(fasta, chrom, start, end, row["strand"], count, row.get("source", "junctions")))
    return result


def write_junctions(path, junctions):
    with open(path, "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "start", "end", "strand", "support", "source", "motif"])
        for j in sorted(junctions):
            writer.writerow([j.chrom, j.start, j.end, j.strand, j.support, j.source, j.motif])


def splice_paths(junctions, lengths, max_nt, max_junctions=3, max_paths=1000, stats=None):
    """Yield bounded local splice paths, including closely spaced microexons.

    Paths combine supported edges, not necessarily read-phased transcripts.
    A separate path cap per first junction bounds combinatorial graph growth.
    """
    stats = stats if stats is not None else Counter()
    grouped = defaultdict(dict)
    for j in junctions:
        key = (j.start, j.end)
        previous = grouped[(j.chrom, j.strand)].get(key)
        if previous is None or previous.support < j.support:
            grouped[(j.chrom, j.strand)][key] = j
        if previous is not None:
            sources = ",".join(sorted(set(previous.source.split(",")) | set(j.source.split(","))))
            grouped[(j.chrom, j.strand)][key] = Junction(j.chrom, j.start, j.end, j.strand, max(previous.support, j.support), sources, j.motif)
    flank = max_nt - 1
    for (chrom, strand), by_coords in sorted(grouped.items()):
        edges = sorted(by_coords.values(), key=lambda j: (j.start, j.end))
        starts = [j.start for j in edges]
        for first in edges:
            stack = [([first], 0)]
            emitted = 0
            while stack:
                if max_paths and emitted >= max_paths:
                    stats["splice_path_roots_capped"] += 1
                    break
                chain, interior_nt = stack.pop()
                segments = [(max(0, chain[0].start - flank), chain[0].start)]
                segments.extend((a.end, b.start) for a, b in zip(chain, chain[1:]))
                segments.append((chain[-1].end, min(lengths[chrom], chain[-1].end + flank)))
                ordered_chain = chain
                if strand == "-":
                    segments.reverse()
                    ordered_chain = list(reversed(chain))
                emitted += 1
                stats["splice_paths_searched"] += 1
                yield chrom, strand, tuple(segments), min(j.support for j in chain), ",".join(j.motif for j in ordered_chain), ",".join(sorted({j.source for j in chain}))
                if len(chain) >= max_junctions:
                    continue
                last = chain[-1]
                low = bisect_left(starts, last.end + 1)
                high = bisect_left(starts, last.end + max_nt - interior_nt - 1)
                for nxt in reversed(edges[low:high]):
                    distance = nxt.start - last.end
                    stack.append((chain + [nxt], interior_nt + distance))
