"""Read transcript exon structures from GTF and GFF3-style Helixer annotations."""
from collections import Counter, defaultdict
import gzip
import re
from urllib.parse import unquote

from .junctions import Junction, motif


_ATTRIBUTE = re.compile(
    r'(?:^|;)\s*([^\s=;]+)\s*(?:=\s*|\s+)'
    r'("(?:[^"\\]|\\.)*"|\'(?:[^\'\\]|\\.)*\'|[^;]*)'
)


def annotation_parents(attributes):
    """Read GTF/GFF2/GFF3 relationships without grouping separate isoforms by gene.

    Split GFF3 Parent lists before URL decoding: an encoded comma is part of an ID.
    Protein-only and bare transcript labels support older predictor output.
    """
    attrs = {}
    for match in _ATTRIBUTE.finditer(attributes):
        key, value = match.groups()
        value = value.strip()
        quoted = len(value) >= 2 and value[0] == value[-1] and value[0] in {'"', "'"}
        if quoted:
            value = re.sub(r'\\([\\"\'])', r'\1', value[1:-1])
        attrs[key.lower()] = value
    for key in ('transcript_id', 'parent', 'protein_id'):
        value = attrs.get(key)
        if value and value != '.':
            parts = value.split(',') if key == 'parent' else [value]
            return [unquote(part.strip()) for part in parts if part.strip()]
    bare = attributes.strip().rstrip(';').strip()
    if re.fullmatch(r'[\w.:%+|/-]+', bare) and bare != '.':
        return [unquote(bare)]
    return []


def read_annotation(path, fasta, lengths):
    """Return exon-adjacency junctions, without pretending annotation is read support.

    Exons are preferred per transcript; CDS intervals are a fallback when exons
    are absent. CDS phase is not trimmed: every junction is searched in all frames,
    preserving codons split across introns. Annotated introns are not BAM-filtered.
    """
    transcripts = defaultdict(lambda: {"exon": set(), "CDS": set()})
    identities = {}
    stats = Counter()
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith('##FASTA'):
                break
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f"annotation line {number}: expected 9 GTF/GFF columns")
            chrom, _, feature, start, end, _, strand, _, attributes = fields
            feature = feature.lower()
            if feature not in {'exon', 'cds'}:
                continue
            feature = 'CDS' if feature == 'cds' else 'exon'
            stats['feature_rows'] += 1
            if chrom not in lengths:
                stats['feature_rows_unknown_contig'] += 1
                continue
            start, end = int(start)-1, int(end)
            if not 0 <= start < end <= lengths[chrom]:
                raise ValueError(f"annotation line {number}: interval outside assembly")
            if strand not in {'+', '-'}:
                raise ValueError(f"annotation line {number}: exon/CDS requires + or - strand")
            parents = annotation_parents(attributes)
            if not parents:
                raise ValueError(f"annotation line {number}: exon/CDS lacks transcript_id or Parent (or protein_id/bare transcript label); gene-only features cannot resolve isoforms")
            for parent in parents:
                identity = (chrom, strand)
                if parent in identities and identities[parent] != identity:
                    raise ValueError(f"annotation transcript {parent!r} appears on conflicting contigs/strands")
                identities[parent] = identity
                transcripts[(chrom, strand, parent)][feature].add((start, end))
    if not transcripts:
        raise ValueError('annotation contains no usable exon/CDS transcripts on assembly contigs')
    junctions = {}
    for (chrom, strand, parent), features in transcripts.items():
        exons = sorted(features['exon'] or features['CDS'])
        stats['transcripts'] += 1
        stats['transcripts_exon' if features['exon'] else 'transcripts_cds_fallback'] += 1
        # Merge touching/overlapping segments; they do not define an intron.
        merged = []
        for a,b in exons:
            if merged and a <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(b, merged[-1][1]))
            else:
                merged.append((a,b))
        for left, right in zip(merged, merged[1:]):
            start, end = left[1], right[0]
            key = (chrom, start, end, strand)
            junctions[key] = Junction(chrom, start, end, strand, 0, 'annotation', motif(fasta, chrom, start, end, strand))
    stats['junctions_retained'] = len(junctions)
    return sorted(junctions.values()), dict(stats)
