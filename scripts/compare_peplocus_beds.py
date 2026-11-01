#!/usr/bin/env python3
"""Independent, strand-aware BED comparison for the carrot integration run."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

import pysam
from Bio.Seq import Seq


def canonical(sequence):
    return sequence.replace('L','I')


def read_bed(path, reverse_negative=False):
    with Path(path).open() as handle:
        for number,line in enumerate(handle,1):
            if not line.strip() or line.startswith(('track','browser','#')):
                continue
            f=line.rstrip('\n').split('\t')
            chrom,start,end,name,strand=f[0],int(f[1]),int(f[2]),f[3],f[5]
            if reverse_negative and strand=='-':name=name[::-1]
            sizes=[int(v) for v in f[10].rstrip(',').split(',')]
            offsets=[int(v) for v in f[11].rstrip(',').split(',')]
            blocks=tuple((start+offset,start+offset+size) for offset,size in zip(offsets,sizes))
            valid=(int(f[9])==len(sizes)==len(offsets) and blocks and strand in {'+','-'}
                   and blocks[0][0]==start and blocks[-1][1]==end
                   and all(a<b for a,b in blocks)
                   and all(a[1]<=b[0] for a,b in zip(blocks,blocks[1:])))
            yield number,name,(canonical(name),chrom,strand,blocks),valid


def checksum(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def compare(root,reference,genome,peptides,prior):
    # This supplied reference has the historical GAPMS negative-label convention.
    # This option is explicit and never auto-detected per record.
    input_classes={canonical(line.strip()) for line in Path(peptides).read_text().splitlines()[1:]}
    index=next((root/'gtf/peplocus/.cache').glob('*.fai'))
    audit=Counter();expected=set();raw=set()
    with pysam.FastaFile(str(genome),filepath_index=str(index)) as fasta, (root/'reference_audit.tsv').open('w') as handle:
        writer=csv.writer(handle,delimiter='\t');writer.writerow(['line','peptide','chrom','strand','blocks','translated','status'])
        for line,name,key,valid in read_bed(reference,reverse_negative=True):
            raw.add(key);_,chrom,strand,blocks=key;observed=''
            valid=valid and chrom in fasta.references and all(0<=a<b<=fasta.get_reference_length(chrom) for a,b in blocks)
            status='invalid_blocks'
            if valid:
                dna=''.join(fasta.fetch(chrom,a,b) for a,b in blocks)
                if strand=='-':dna=str(Seq(dna).reverse_complement())
                if len(dna)!=3*len(name):status='coding_length_mismatch'
                else:
                    observed=str(Seq(dna).translate())
                    status='exact' if observed==name else 'il_equivalent' if canonical(observed)==canonical(name) else 'sequence_mismatch'
            audit[status]+=1
            if status in {'exact','il_equivalent'} and key[0] in input_classes:expected.add(key)
            writer.writerow([line,name,chrom,strand,json.dumps(blocks),observed,status])
    summary={'reference':str(reference),'reference_sha256':checksum(reference),
             'reference_label_convention':'negative-strand names reversed back to N-to-C',
             'reference_audit':dict(audit),'reference_evaluable_loci':len(expected),
             'reference_evaluable_classes':len({k[0] for k in expected}),
             'input_classes':len(input_classes),'modes':{}}
    for mode in ['gtf','bam','workflow']:
        bed=root/mode/'peptides_mapped.bed';igv=root/mode/'peptides_mapped.igv.bed'
        found={key for _,_,key,valid in read_bed(bed) if valid}
        display={key for _,_,key,valid in read_bed(igv,reverse_negative=True) if valid}
        if found!=display:raise ValueError(f'{mode}: IGV differs from standard placements')
        counts=Counter(k[0] for k in found);correct=expected&found;missing=expected-found
        report={'reported_loci':len(found),'mapped_classes':len(counts),
                'unmapped_classes':len(input_classes-counts.keys()),
                'unique_classes':sum(n==1 for n in counts.values()),
                'multiple_classes':sum(n>1 for n in counts.values()),
                'correct_reference_loci':len(correct),'missed_reference_loci':len(missing),
                'reference_locus_recall':len(correct)/len(expected),
                'reference_peptides_recovered':len({k[0] for k in correct}),
                'additional_loci':len(found-expected),'raw_reference_coordinate_overlap':len(found&raw),
                'igv_round_trip_equal':True,
                'engine_summary':json.loads((root/mode/'peplocus/summary.json').read_text())}
        if mode=='gtf':
            old={key for _,_,key,valid in read_bed(prior) if valid}
            report['previous_standalone_peplocus_loci']=len(old)
            report['identical_to_previous_standalone_peplocus']=found==old
        if mode=='workflow':
            direct={key for _,_,key,valid in read_bed(root/'gtf/peptides_mapped.bed') if valid}
            report['identical_to_standalone_gapms_map']=found==direct
        summary['modes'][mode]=report
        with (root/f'{mode}_missed_reference_loci.tsv').open('w') as handle:
            writer=csv.writer(handle,delimiter='\t');writer.writerow(['peptide_class','chrom','strand','blocks'])
            for peptide,chrom,strand,blocks in sorted(missing):writer.writerow([peptide,chrom,strand,json.dumps(blocks)])
    (root/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines=['# Carrot GAPMS/PepLocus integration comparison','','Reference: `'+str(reference)+'`. Negative-strand display labels were restored before comparison. Locus matching uses I/L-equivalent peptide, chromosome, strand, and complete BED blocks.','','| Mode | Mapped queries | Loci | Recovered reference loci | Missed | Recall | Multiple queries | Seconds |','|---|---:|---:|---:|---:|---:|---:|---:|']
    for mode,m in summary['modes'].items():
        lines.append(f"| {mode} | {m['mapped_classes']:,} | {m['reported_loci']:,} | {m['correct_reference_loci']:,} | {m['missed_reference_loci']:,} | {m['reference_locus_recall']:.4%} | {m['multiple_classes']:,} | {m['engine_summary']['elapsed_seconds']:.3f} |")
    lines+=['',f"Reference audit: {dict(audit)}. The recall denominator is {len(expected):,} sequence-consistent reference loci, not every original BED row.",'','Additional sequence-compatible loci are not automatically false positives; the protein-derived reference does not exhaust genomic placements. All IGV outputs reproduce standard loci and query identities when their display convention is undone. GTF output is checked against the prior standalone PepLocus run, and workflow output against gapms map. Full configuration and equality results are in comparison.json. BAM timing includes extraction; GTF/workflow mapping timing excludes protein scoring.','']
    (root/'REPORT.md').write_text('\n'.join(lines))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('root',type=Path)
    parser.add_argument('--reference',required=True,type=Path);parser.add_argument('--genome',required=True,type=Path)
    parser.add_argument('--peptides',required=True,type=Path);parser.add_argument('--prior-peplocus',required=True,type=Path)
    args=parser.parse_args();compare(args.root,args.reference,args.genome,args.peptides,args.prior_peplocus)
