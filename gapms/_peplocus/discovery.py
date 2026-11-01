"""Optional bounded, annotation-free discovery of single-intron peptide mappings."""
from array import array
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
import logging

import ahocorasick

from .sequence import canonical, reverse_complement, translate


LOG = logging.getLogger("peplocus")


def discover(fasta, lengths, peptides, peptide_ids, table=1, il_equivalent=True, anchor=4, min_intron=20, max_intron=500000, max_seed_hits=2000, max_pairs=10000, chunk_size=2000000, stats=None, limited=None, cache_mb=256):
    """Join exact terminal amino-acid seeds and test all nucleotide splice phases.

    Exhaustive only within the configured single-intron, canonical-motif, anchor,
    intron-length, and repeat/candidate limits. No protein prediction is involved.
    """
    stats = stats if stats is not None else Counter()
    limited = limited if limited is not None else set()
    targets = [pid for pid in peptide_ids if len(peptides[pid]) >= 2 * anchor]
    stats["novel_eligible_peptides"] = len(targets)
    if not targets:
        return
    seeds = {}
    terminal = {}
    for pid in targets:
        seq = canonical(peptides[pid], il_equivalent)
        left, right = seq[:anchor], seq[-anchor:]
        seeds.setdefault(left, len(seeds))
        seeds.setdefault(right, len(seeds))
        terminal[pid] = (seeds[left], seeds[right])
    automaton = ahocorasick.Automaton()
    for word, sid in seeds.items():
        automaton.add_word(word, sid)
    automaton.make_automaton()
    for chrom, length in lengths.items():
        LOG.info("Novel discovery: collecting seed positions on %s (%s bp)", chrom, length)
        # Small random faidx reads are expensive on shared filesystems. Cache
        # a contig only when the forward/reverse/transient copies fit the budget;
        # larger chromosomes retain indexed on-demand access.
        cached_forward = fasta.fetch(chrom).upper() if length * 3 <= cache_mb * 1024 * 1024 else None
        if cached_forward is not None:
            stats["novel_contigs_cached"] += 1
        positions = defaultdict(lambda: array("Q"))
        capped = set()
        overlap = 3 * anchor - 1
        for start in range(0, length, chunk_size):
            core_end = min(length, start + chunk_size)
            end = min(length, core_end + overlap)
            forward = cached_forward[start:end] if cached_forward is not None else fasta.fetch(chrom, start, end).upper()
            for strand, dna in (("+", forward), ("-", reverse_complement(forward))):
                for frame in range(3):
                    protein = canonical(translate(dna[frame:], table), il_equivalent)
                    for last, sid in automaton.iter(protein):
                        local = frame + 3 * (last - anchor + 1)
                        genomic = start + local if strand == "+" else end - local - 3 * anchor
                        if not start <= genomic < core_end:
                            continue
                        key = (strand, sid)
                        if max_seed_hits and len(positions[key]) >= max_seed_hits:
                            capped.add(key)
                            continue
                        position = genomic if strand == "+" else length - genomic - 3 * anchor
                        positions[key].append(position)
        stats["novel_seed_keys_capped"] += len(capped)
        LOG.info("Novel discovery: %s retained seed positions on %s; %s capped keys", sum(map(len, positions.values())), chrom, len(capped))
        for key, values in positions.items():
            positions[key] = array("Q", sorted(values))
        for strand in ("+", "-"):
            oriented = cached_forward if strand == "+" else reverse_complement(cached_forward) if cached_forward is not None else None
            def fetch(a, b):
                if oriented is not None:
                    return oriented[a:b]
                if strand == "+":
                    return fasta.fetch(chrom, a, b).upper()
                return reverse_complement(fasta.fetch(chrom, length - b, length - a).upper())
            for pid in targets:
                left_seed, right_seed = terminal[pid]
                if (strand, left_seed) in capped or (strand, right_seed) in capped:
                    limited.add(pid)
                lefts = positions.get((strand, left_seed), ())
                rights = positions.get((strand, right_seed), ())
                peptide = peptides[pid]
                nt_length = len(peptide) * 3
                expected = canonical(peptide, il_equivalent)
                # Verify each seed's surrounding peptide context before joining
                # loci. Most random 4-aa anchors cannot extend to a valid splice
                # boundary. This avoids the Cartesian product of such anchors.
                right_groups = defaultdict(list)
                right_windows = {}
                for right in rights:
                    end = right + 3 * anchor
                    window_start = end - nt_length
                    if window_start < 0:
                        continue
                    window = fetch(window_start, end)
                    protein = canonical(translate(window, table), il_equivalent)
                    stats["novel_seed_contexts_checked"] += 1
                    shared = 0
                    for observed, wanted in zip(reversed(protein), reversed(expected)):
                        if observed != wanted:
                            break
                        shared += 1
                    minimum_cut = max(3 * anchor, 3 * (len(peptide) - shared) - 2)
                    accepted = False
                    for cut in range(minimum_cut, nt_length - 3 * anchor + 1):
                        acceptor = window[cut-2:cut]
                        if acceptor in {"AG", "AC"}:
                            right_groups[(cut, acceptor)].append(end)
                            accepted = True
                    if accepted:
                        right_windows[end] = window
                pairs = 0
                exhausted = False
                for start in lefts:
                    if start + nt_length > length:
                        continue
                    left_window = fetch(start, start + nt_length)
                    protein = canonical(translate(left_window, table), il_equivalent)
                    stats["novel_seed_contexts_checked"] += 1
                    shared = 0
                    for observed, wanted in zip(protein, expected):
                        if observed != wanted:
                            break
                        shared += 1
                    maximum_cut = min(nt_length - 3 * anchor, 3 * shared + 2)
                    for cut in range(3 * anchor, maximum_cut + 1):
                        donor_motif = left_window[cut:cut+2]
                        if donor_motif not in {"GT", "GC", "AT"}:
                            continue
                        acceptor_motif = "AC" if donor_motif == "AT" else "AG"
                        ends = right_groups.get((cut, acceptor_motif), ())
                        low = start + nt_length + min_intron
                        high = start + nt_length + max_intron
                        lo, hi = bisect_left(ends, low), bisect_right(ends, high)
                        for end in ends[lo:hi]:
                            if max_pairs and pairs >= max_pairs:
                                limited.add(pid)
                                stats["novel_pair_groups_capped"] += 1
                                exhausted = True
                                break
                            pairs += 1
                            stats["novel_anchor_pairs_tested"] += 1
                            donor = start + cut
                            acceptor = end - (nt_length - cut)
                            right_window = right_windows[end]
                            matched = translate(left_window[:cut] + right_window[cut:], table)
                            if canonical(matched, il_equivalent) != expected:
                                continue
                            if strand == "+":
                                blocks = ((start, donor), (acceptor, end))
                            else:
                                blocks = ((length - end, length - acceptor), (length - donor, length - start))
                            m = donor_motif + "-" + acceptor_motif
                            stats["novel_hits_before_deduplication"] += 1
                            yield (pid, chrom, strand, blocks, matched, "novel", 0, m)
                        if exhausted:
                            break
                    if exhausted:
                        break
        oriented = None
        cached_forward = None
        LOG.info("Novel discovery: finished %s", chrom)
