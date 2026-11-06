"""Bounded structural support from optional RNA/protein evidence."""

from collections import defaultdict
import json, math, bisect
from pathlib import Path
from .intervals import merge, intersection_length


class Support:
    def __init__(self, rna=None, proteins=(), min_rna_templates=3):
        self.min_rna_templates = min_rna_templates
        self.rna = {}
        self.rna_blocks = {}
        self.rna_internal = {}
        self.rna_exon_counts = {}
        self.rna_exon_counts_all = {}
        self.rna_chain_counts = {}
        self.rna_block_starts = {}
        self.rna_block_ends = {}
        self.protein_introns = set()
        self.protein_exon_frames = set()
        self.protein_blocks = {}
        self.rna_bins = defaultdict(set)
        self.protein_bins = defaultdict(set)
        if rna:
            data = (
                json.loads(Path(rna).read_text())
                if isinstance(rna, (str, Path))
                else rna
            )
            if data.get("coordinate_system") != "0-based-half-open":
                raise ValueError(
                    "RNA junctions must specify 0-based-half-open coordinates"
                )
            rb=defaultdict(list);ri=defaultdict(list)
            for field,target in [('alignment_blocks',rb),('observed_internal_exons',ri)]:
                for chrom,a,b,strand in data.get(field,[]):
                    if chrom not in data['references'] or strand not in ('+','-') or not isinstance(a,int) or not isinstance(b,int) or not 0<=a<b<=data['references'][chrom]:raise ValueError('Invalid RNA observed-block geometry')
                    target[chrom,strand].append((a,b))
            self.rna_blocks={key:merge(v) for key,v in rb.items()}
            self.rna_internal={key:tuple(sorted(set(v))) for key,v in ri.items()}
            self.rna_block_starts={k:tuple(a for a,b in v) for k,v in self.rna_blocks.items()}
            self.rna_block_ends={k:tuple(b for a,b in v) for k,v in self.rna_blocks.items()}
            for field,target in [('observed_internal_exon_counts',self.rna_exon_counts),('rna_chain_junctions',self.rna_chain_counts),('observed_internal_exon_counts_all',self.rna_exon_counts_all)]:
                for chrom,a,b,strand,n in data.get(field,[]):
                    if chrom not in data['references'] or strand not in ('+','-') or not 0<=a<b<=data['references'][chrom] or not isinstance(n,int) or n<1:
                        raise ValueError('Invalid RNA structure count')
                    key=(chrom,strand,a,b)
                    target[key]=max(target.get(key,0),n)
            if 'observed_internal_exon_counts_all' in data and any(self.rna_exon_counts_all.get(k,0)<n for k,n in self.rna_exon_counts.items()):
                raise ValueError('Unfiltered RNA exon counts must include every filtered exon count')
            for chrom, a, b, strand, n in data["junctions"]:
                if (
                    strand not in ("+", "-", ".")
                    or not 0 <= a < b
                    or not isinstance(n, int)
                    or n < 1
                ):
                    raise ValueError("Invalid RNA junction record")
                key = (chrom, strand, a, b)
                # Duplicate rows are not independent observations.
                self.rna[key] = max(self.rna.get(key, 0), n)
                self.rna_bins[chrom, a // 10000].add(key)
        blocks = defaultdict(list)
        seen = set()
        for model in proteins:
            if model.key in seen:
                continue
            seen.add(model.key)
            used=-model.phase
            for iv in (model.cds if model.strand=='+' else reversed(model.cds)):
                self.protein_exon_frames.add((model.contig,model.strand,iv,used%3))
                used+=iv[1]-iv[0]

            blocks[model.contig, model.strand].extend(model.cds)
            self.protein_introns.update(
                (model.contig, model.strand, a, b) for a, b in model.introns
            )
        for chrom, strand, a, b in self.protein_introns:
            self.protein_bins[chrom, strand, a // 10000].add((a, b))
        self.protein_blocks = {key: merge(v) for key, v in blocks.items()}
        self.protein_starts = {
            k: [a for a, b in v] for k, v in self.protein_blocks.items()
        }
        self.protein_ends = {
            k: [b for a, b in v] for k, v in self.protein_blocks.items()
        }

    def rna_coding_overlap(self, model):
        key = (model.contig, model.strand)
        blocks = self.rna_blocks.get(key, ())
        if not blocks:
            return 0
        left = bisect.bisect_right(self.rna_block_ends[key], model.cds[0][0])
        right = bisect.bisect_left(self.rna_block_starts[key], model.cds[-1][1])
        return intersection_length(model.cds, blocks[left:right])

    def rna_count(self, model, iv):
        a, b = iv
        return max(
            self.rna.get((model.contig, model.strand, a, b), 0),
            self.rna.get((model.contig, ".", a, b), 0),
        )

    def rna_score(self, model, introns=None):
        ivs = model.introns if introns is None else tuple(introns)
        if not ivs:
            return 0.0
        return sum(
            (n / (n + 3.0) if n >= self.min_rna_templates else 0.0)
            for n in (self.rna_count(model, iv) for iv in ivs)
        ) / len(ivs)

    def protein_score(self, model, regions=None):
        intervals = model.cds if regions is None else tuple(regions)
        size = sum(b - a for a, b in intervals)
        key = (model.contig, model.strand)
        blocks = self.protein_blocks.get(key, ())
        if intervals and blocks:
            first = bisect.bisect_right(
                self.protein_ends[key], min(a for a, b in intervals)
            )
            last = bisect.bisect_left(
                self.protein_starts[key], max(b for a, b in intervals)
            )
            blocks = blocks[first:last]
        coverage = intersection_length(merge(intervals), blocks) / size if size else 0.0
        ivs = model.introns
        junctions = (
            sum(
                (model.contig, model.strand, a, b) in self.protein_introns
                for a, b in ivs
            )
            / len(ivs)
            if ivs
            else coverage
        )
        return (coverage + junctions) / 2

    def supported_junction(self, model, iv, predictor_votes=None):
        a, b = iv
        return (
            (predictor_votes or {}).get((model.contig, model.strand, a, b), 0) >= 2
            or self.rna_count(model, iv) >= self.min_rna_templates
            or (model.contig, model.strand, a, b) in self.protein_introns
        )

    def graph_junctions(self, models):
        if not models:
            return []
        m = models[0]
        lo = min(x.cds[0][0] for x in models)
        hi = max(x.cds[-1][1] for x in models)
        keys = set()
        out = []
        for k in range(lo // 10000, (hi - 1) // 10000 + 1):
            keys.update(self.rna_bins[m.contig, k])
        for chrom, strand, a, b in sorted(keys):
            n = self.rna[chrom, strand, a, b]
            if (
                strand in (".", m.strand)
                and lo <= a < b <= hi
                and n >= self.min_rna_templates
            ):
                out.append((a, b, n / (n + 3.0), "rna"))
        keys = set()
        for k in range(lo // 10000, (hi - 1) // 10000 + 1):
            keys.update(self.protein_bins[m.contig, m.strand, k])
        out.extend((a, b, 0.5, "protein") for a, b in sorted(keys) if lo <= a < b <= hi)
        return out

    def graph_segments(self, models):
        if not models:
            return []
        m = models[0]
        lo = min(x.cds[0][0] for x in models)
        hi = max(x.cds[-1][1] for x in models)
        key = (m.contig, m.strand)
        blocks = self.protein_blocks.get(key, ())
        rna_segments=[(max(lo,a),min(hi,b),'rna') for a,b in self.rna_internal.get(key,()) if max(lo,a)<min(hi,b)]
        if not blocks:
            return rna_segments
        first = bisect.bisect_right(self.protein_ends[key], lo)
        last = bisect.bisect_left(self.protein_starts[key], hi)
        return rna_segments + [
            (max(lo, a), min(hi, b), "protein")
            for a, b in blocks[first:last]
            if max(lo, a) < min(hi, b)
        ]
