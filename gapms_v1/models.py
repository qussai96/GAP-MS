"""Source-aware coding paths in zero-based, half-open genomic coordinates."""
from dataclasses import dataclass, field
import gzip
import json
from pathlib import Path
import re
from urllib.parse import unquote

from gapms._peplocus.annotation import _ATTRIBUTE, annotation_parents
from gapms._peplocus.sequence import join_blocks, translate, project
from .common import stable_id


def attributes(text):
    result = {}
    for match in _ATTRIBUTE.finditer(text):
        key, value = match.groups()
        result[key.lower()] = unquote(value.strip().strip('"\''))
    return result


@dataclass
class Model:
    id: str
    source: str
    original_id: str
    gene: str
    chrom: str
    strand: str
    exons: list = field(default_factory=list)
    cds: list = field(default_factory=list)
    phases: dict = field(default_factory=dict)
    protein: str = ''
    coding: list = field(default_factory=list)
    issues: list = field(default_factory=list)
    complete_start: bool = False
    complete_stop: bool = False

    @property
    def intervals(self):
        return self.exons or self.cds

    @property
    def junctions(self):
        return tuple((a[1], b[0]) for a, b in zip(self.cds, self.cds[1:]) if a[1] < b[0])

    def project(self, start, end):
        return project(start, end, self.coding, self.strand)

    def signature(self):
        return (self.chrom, self.strand, tuple(self.coding))


def load_models(path, source, reference, lengths, genetic_code=1, allow_unstranded=False, namespace_loci=False, skipped_features=None):
    """Preserve isoforms and GFF Parent lists; trim only the first CDS phase."""
    opener = gzip.open if str(path).endswith('.gz') else open
    models = {}
    gene_parents = {}
    rows = []
    with opener(path, 'rt') as stream:
        for number, line in enumerate(stream, 1):
            if line.startswith('##FASTA'):
                break
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 9:
                raise ValueError(f'{path}:{number}: expected 9 annotation columns')
            chrom, _, feature, left, right, _, strand, phase, text = parts
            feature = feature.lower()
            attrs = attributes(text)
            if feature in {'mrna', 'transcript'}:
                tid = attrs.get('id') or attrs.get('transcript_id')
                if tid:
                    gene_parents[tid] = attrs.get('gene_id') or attrs.get('parent') or tid
            if feature not in {'exon', 'cds'}:
                continue
            left, right = int(left) - 1, int(right)
            if chrom not in lengths or not 0 <= left < right <= lengths[chrom]:
                raise ValueError(f'{path}:{number}: feature outside assembly or unknown contig {chrom}')
            if strand not in {'+', '-'} and not (allow_unstranded and strand=='.' and feature=='exon'):
                raise ValueError(f'{path}:{number}: exon/CDS strand must be + or -')
            parents = annotation_parents(text)
            if not parents and namespace_loci and feature=='exon':
                if skipped_features is not None:skipped_features.append(dict(path=str(path),line=number,feature=feature,reason='orphan_exon_without_transcript_parent_not_a_coding_path'))
                continue
            if not parents:
                raise ValueError(f'{path}:{number}: missing transcript identity')
            if phase not in {'.', '0', '1', '2'}:
                raise ValueError(f'{path}:{number}: invalid CDS phase {phase}')
            rows.append((parents, chrom, strand, feature, left, right, phase, attrs))
    for parents, chrom, strand, feature, left, right, phase, attrs in rows:
        for tid in parents:
            model = models.setdefault((chrom,strand,tid) if namespace_loci else tid, Model(stable_id('M', source, (chrom,strand,tid) if namespace_loci else tid), source, tid,
                    attrs.get('gene_id') or gene_parents.get(tid, tid), chrom, strand))
            if (model.chrom, model.strand) != (chrom, strand):
                raise ValueError(f'{path}: transcript {tid} has conflicting contigs or strands')
            interval = (left, right)
            if feature == 'exon':
                model.exons.append(interval)
            else:
                if interval in model.phases and model.phases[interval] != phase:
                    model.issues.append('conflicting_duplicate_CDS_phase')
                model.cds.append(interval)
                model.phases[interval] = phase
    if not models:
        raise ValueError(f'{path}: no exon/CDS transcript models')
    for model in models.values():
        model.exons = sorted(set(model.exons))
        model.cds = sorted(set(model.cds))
        for intervals in [model.exons, model.cds]:
            if any(a[1] > b[0] for a, b in zip(intervals, intervals[1:])):
                raise ValueError(f'{path}: overlapping features within transcript {model.original_id}')
        if not model.cds:
            model.issues.append('no_CDS_annotation')
            continue
        coding = list(model.cds if model.strand == '+' else reversed(model.cds))
        first_phase = model.phases[coding[0]]
        if first_phase == '.':
            model.issues.append('first_CDS_phase_unknown_assumed_zero')
        trim = int(first_phase) if first_phase != '.' else 0
        remaining = trim
        trimmed = []
        for left, right in coding:
            removed = min(remaining, right-left)
            remaining -= removed
            if model.strand == '+':left += removed
            else:right -= removed
            if left < right:trimmed.append((left, right))
        model.coding = trimmed
        cumulative = coding[0][1] - coding[0][0] - trim
        for interval in coding[1:]:
            p = model.phases[interval]
            if p != '.' and int(p) != (-cumulative) % 3:
                model.issues.append('inconsistent_internal_CDS_phase')
            cumulative += interval[1] - interval[0]
        dna = join_blocks(reference, model.chrom, model.strand, model.coding)
        aa = translate(dna, genetic_code)
        model.complete_start = dna[:3].upper() == 'ATG'
        model.complete_stop = len(dna) % 3 == 0 and aa.endswith('*')
        model.protein = aa[:-1] if aa.endswith('*') else aa
        if '*' in model.protein:model.issues.append('internal_stop')
        if len(dna) % 3:model.issues.append('partial_terminal_codon')
    return sorted(models.values(), key=lambda m: m.id)


def write_models(models, path):
    """Write CDS paths with explicit source and stable model IDs; no input edits."""
    with Path(path).open('w') as stream:
        for model in sorted(models, key=lambda m: (m.chrom, min(m.intervals), m.strand, m.id)):
            intervals = model.intervals
            left = min(a for a,b in intervals); right = max(b for a,b in intervals)
            # JSON string escaping is compatible with quoted GTF attributes.
            attrs = '; '.join(k+' '+json.dumps(str(v)) for k,v in [
                ('gene_id', model.source+':'+model.gene), ('transcript_id', model.id),
                ('original_transcript_id',model.original_id),('annotation_source',model.source)])+';'
            stream.write(f'{model.chrom}\tGAPMS_v1\ttranscript\t{left+1}\t{right}\t.\t{model.strand}\t.\t{attrs}\n')
            for a,b in model.exons:
                stream.write(f'{model.chrom}\tGAPMS_v1\texon\t{a+1}\t{b}\t.\t{model.strand}\t.\t{attrs}\n')
            for a,b in model.cds:
                phase = model.phases.get((a,b), '.')
                stream.write(f'{model.chrom}\tGAPMS_v1\tCDS\t{a+1}\t{b}\t.\t{model.strand}\t{phase}\t{attrs}\n')


def write_proteins(models, path):
    with Path(path).open('w') as stream:
        for model in sorted(models, key=lambda m:m.id):
            if model.protein:
                stream.write(f'>{model.id}\n{model.protein}\n')
