"""Bounded, strand-aware exon-segment / splice-junction DAG reconstruction.

Existing coding coverage and explicit optional RNA/protein alignment segments
can form exon segments. RNA segments require observed internal-exon geometry.
Path search is bounded and reports pruning. Final ORFs are checked on the genome.
"""

from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass
from .models import Transcript


@dataclass(frozen=True)
class Edge:
    start: int
    end: int
    kind: str
    sources: tuple
    strength: float


@dataclass(frozen=True)
class PathState:
    score: float
    blocks: tuple
    bases: int
    junctions: tuple
    start: int
    observed_peptides: tuple = ()


def loci(models):
    """Deterministic envelope components, separated by contig and strand."""
    grouped = defaultdict(list)
    for m in models:
        grouped[m.contig, m.strand].append(m)
    for key, mm in sorted(grouped.items()):
        component = []
        right = -1
        for m in sorted(mm, key=lambda m: (m.cds[0][0], m.cds[-1][1], m.uid)):
            if component and m.cds[0][0] >= right:
                yield component
                component = []
                right = -1
            component.append(m)
            right = max(right, m.cds[-1][1])
        if component:
            yield component


def build_graph(models, extra_junctions=(), extra_segments=()):
    if not models:
        return {}, [], []
    chrom, strand = models[0].contig, models[0].strand
    if any((m.contig, m.strand) != (chrom, strand) for m in models):
        raise ValueError("A graph cannot cross contigs or strands")
    bounds = {p for m in models for iv in m.cds for p in iv}
    coverage = defaultdict(list)
    junctions = defaultdict(set)
    for m in models:
        for a, b in m.cds:
            coverage[a].append((1, m.source))
            coverage[b].append((-1, m.source))
        for iv in m.introns:
            junctions[iv].add(m.source)
    extra = {}
    lo = min(bounds)
    hi = max(bounds)
    for a, b, label in extra_segments:
        a, b = max(lo, a), min(hi, b)
        if a < b:
            bounds.update((a, b))
            coverage[a].append((1, label))
            coverage[b].append((-1, label))
    for a, b, strength, label in extra_junctions:
        if lo <= a < b <= hi:
            bounds.update((a, b))
            extra[a, b] = max(extra.get((a, b), 0.0), float(strength))
            junctions[a, b].add(label)
    bounds = sorted(bounds)
    active = defaultdict(int)
    graph = defaultdict(list)
    direction = 1 if strand == "+" else -1
    for a, b in zip(bounds, bounds[1:]):
        for delta, source in coverage[a]:
            active[source] += delta
        ss = tuple(sorted(s for s, n in active.items() if n > 0))
        if ss:
            start, end = (a, b) if direction == 1 else (b, a)
            graph[start].append(Edge(start, end, "coding", ss, len(ss)))
    for (a, b), ss in sorted(junctions.items()):
        start, end = (a, b) if direction == 1 else (b, a)
        graph[start].append(
            Edge(
                start,
                end,
                "splice",
                tuple(sorted(ss)),
                len(ss) + extra.get((a, b), 0.0),
            )
        )
    for edges in graph.values():
        edges.sort(key=lambda e: (direction * e.end, e.kind, e.sources))
    starts = sorted(
        {m.cds[0][0] if direction == 1 else m.cds[-1][1] for m in models},
        reverse=direction == -1,
    )
    ends = sorted(
        {m.cds[-1][1] if direction == 1 else m.cds[0][0] for m in models},
        reverse=direction == -1,
    )
    return dict(graph), starts, ends


def _extend(state, edge, strand, molecular_corrections=None):
    a,b=sorted((edge.start,edge.end))
    correction=(molecular_corrections or {}).get((edge.kind,a,b,state.bases%3),0.)
    if edge.kind == "splice":
        if not state.blocks:
            return None
        terminal = state.blocks[-1][1] if strand == "+" else state.blocks[-1][0]
        if terminal != edge.start:
            return None
        return PathState(
            state.score + 0.2 * (edge.strength + correction),
            state.blocks,
            state.bases,
            state.junctions + ((min(edge.start, edge.end), max(edge.start, edge.end)),),
            state.start,
            state.observed_peptides,
        )
    a, b = sorted((edge.start, edge.end))
    blocks = list(state.blocks)
    # An exon cannot silently traverse a junction: only contiguous coding edges merge.
    if blocks and (
        (strand == "+" and blocks[-1][1] == a) or (strand == "-" and blocks[-1][0] == b)
    ):
        blocks[-1] = (min(a, blocks[-1][0]), max(b, blocks[-1][1]))
    else:
        blocks.append((a, b))
    return PathState(
        state.score + (b - a) / 1000.0 * (edge.strength + correction),
        tuple(blocks),
        state.bases + b - a,
        state.junctions,
        state.start,
        state.observed_peptides,
    )


def reconstruct(
    models,
    genome,
    extra_junctions=(),
    beam=8,
    max_candidates=32,
    max_nodes=2000,
    max_starts=32,
    extra_segments=(),
    molecular_corrections=None,
    peptide_observations=(),
):
    if min(beam, max_candidates, max_nodes, max_starts) < 1:
        raise ValueError("Graph bounds must be positive")
    if not models:
        return [], dict(nodes=0, edges=0, paths_pruned=0, candidates=0)
    graph, starts, ends = build_graph(models, extra_junctions, extra_segments)
    strand = models[0].strand
    chrom = models[0].contig
    # A path with matching geometry/frame uses exactly these peptide bases.
    # Validate their translation once, rather than retranslating whole partial
    # ORFs and reassessing all preceding introns at each search extension.
    if peptide_observations:
        from Bio.Seq import Seq
        validated=[]
        for peptide in peptide_observations:
            if (peptide.contig,peptide.strand)!=(chrom,strand):continue
            if peptide.blocks[0][0]<min(m.cds[0][0] for m in models) or peptide.blocks[-1][1]>max(m.cds[-1][1] for m in models):continue
            dna=''.join(genome.fetch(peptide.contig,a,b) for a,b in peptide.blocks)
            if peptide.strand=='-':dna=str(Seq(dna).reverse_complement())
            if str(Seq(dna).translate()).replace('I','L')==peptide.sequence:
                validated.append(peptide)
        peptide_observations=tuple(validated)
    direction = 1 if strand == "+" else -1
    nodes = sorted(
        {p for v in graph.values() for e in v for p in (e.start, e.end)},
        reverse=direction == -1,
    )
    report = dict(
        nodes=len(nodes),
        edges=sum(map(len, graph.values())),
        paths_pruned=0,
        candidates=0,
        skipped=None,
    )
    if len(nodes) > max_nodes:
        report["skipped"] = "node_limit"
        return [], report
    if len(starts) > max_starts:
        report["skipped"] = "start_limit"
        return [], report
    states = defaultdict(list)
    endset = set(ends)
    known = {m.key for m in models}
    found = {}
    for start in starts:
        states[start].append(PathState(0.0, (), 0, (), start))
    for pos in nodes:
        # Retain a bounded beam separately for each initial boundary and codon phase.
        groups = defaultdict(list)
        for state in states.pop(pos, []):
            groups[state.start, state.bases % 3].append(state)
        current = []
        for key, ss in sorted(groups.items()):
            unique = {(s.blocks, s.junctions): s for s in ss}
            ordered = sorted(
                unique.values(), key=lambda s: (-s.score, s.blocks, s.junctions)
            )
            report["paths_pruned"] += max(0, len(ordered) - beam)
            current.extend(ordered[:beam])
        for state in current:
            if pos in endset and state.blocks and state.bases % 3 == 0:
                # End only after coding, never after a splice into a terminal boundary.
                blocks = tuple(sorted(state.blocks))
                terminal = blocks[-1][1] if direction == 1 else blocks[0][0]
                if terminal == pos:
                    m = Transcript(
                        chrom, strand, blocks, "graph", "graph", 0, blocks, "", "graph"
                    )
                    if m.key not in known and m.key not in found:
                        complete, canonical, _ = genome.assess(m)
                        if complete and canonical:
                            found[m.key] = (state.score, m)
            for edge in graph.get(pos, ()):
                if edge.kind == "splice":
                    a, b = sorted((edge.start, edge.end))
                    # Test splice grammar on the actual reference strand via a tiny model.
                    if b - a < 4:
                        continue
                    from Bio.Seq import Seq

                    left, right = genome.fetch(chrom, a, a + 2), genome.fetch(
                        chrom, b - 2, b
                    )
                    if strand == "-":
                        left, right = str(Seq(right).reverse_complement()), str(
                            Seq(left).reverse_complement()
                        )
                    if (left, right) not in {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}:
                        continue
                for state in current:
                    new = _extend(state, edge, strand, molecular_corrections)
                    if new is not None:
                        if peptide_observations and edge.kind=='coding':
                            # Reward only a complete, frame/translation-compatible
                            # peptide chain. Its component junctions are not free
                            # independent observations that can be recombined.
                            from .peptides import peptide_offset
                            observed=set(new.observed_peptides)
                            relevant=[p for p in peptide_observations if p.sequence not in observed
                                and new.blocks and min(a for a,b in new.blocks)<=p.blocks[0][0]
                                and max(b for a,b in new.blocks)>=p.blocks[-1][1]]
                            if relevant:
                                partial=Transcript(chrom,strand,tuple(sorted(new.blocks)),'graph','partial',0)
                                bonus=0.
                                for peptide in relevant:
                                    if peptide_offset(partial,peptide) is not None:
                                        observed.add(peptide.sequence)
                                        bonus+=3.*len(peptide.sequence)/1000.+.2*(len(peptide.blocks)-1)
                                if bonus:
                                    new=PathState(new.score+bonus,new.blocks,new.bases,new.junctions,
                                                  new.start,tuple(sorted(observed)))
                        states[edge.end].append(new)
    ordered = sorted(found.values(), key=lambda x: (-x[0], x[1].key))
    report["candidates_before_limit"] = len(ordered)
    report["candidates"] = min(len(ordered), max_candidates)
    report["paths_pruned"] += max(0, len(ordered) - max_candidates)
    return [m for score, m in ordered[:max_candidates]], report
