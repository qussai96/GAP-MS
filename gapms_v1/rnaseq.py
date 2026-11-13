"""RNA transcript reconstruction with read-only assembly access."""
from pathlib import Path
import os
import shutil
import json
from collections import defaultdict

import pysam
from peplocus.sequence import join_blocks
from .common import run_command, stable_id, write_table
from .models import load_models


def validate_bam(path, lengths):
    with pysam.AlignmentFile(str(path), 'rb') as bam:
        if bam.header.to_dict().get('HD', {}).get('SO') != 'coordinate':
            raise ValueError('RNA reconstruction requires a coordinate-sorted BAM (HD:SO:coordinate)')
        if not set(bam.references) & set(lengths):
            raise ValueError('BAM and assembly have no matching contig names')
        for name, length in zip(bam.references, bam.lengths):
            if name in lengths and lengths[name] != length:
                raise ValueError(f'BAM/assembly length mismatch for {name}')
    # BAM is read sequentially by junction extraction and StringTie; no input index is created.



def canonicalize_transcripts(models, path, provenance):
    """Remove StringTie scheduling-dependent names/order, preserving gene groups."""
    groups=defaultdict(set)
    for m in models:groups[m.chrom,m.strand,m.gene].add(tuple(m.intervals))
    genes={key:stable_id('RG',key[0],key[1],sorted(paths)) for key,paths in groups.items()}
    unique={};records=[]
    for m in models:
        old_transcript=m.original_id;old_gene=m.gene
        gene=genes[m.chrom,m.strand,old_gene]
        transcript=stable_id('RT',gene,m.chrom,m.strand,m.intervals)
        records.append(dict(original_transcript_id=old_transcript,original_gene_id=old_gene,
                            transcript_id=transcript,gene_id=gene,chrom=m.chrom,strand=m.strand,exon_blocks=m.intervals))
        m.original_id=transcript;m.gene=gene;m.id=stable_id('M',m.source,transcript)
        unique.setdefault(transcript,m)
    result=sorted(unique.values(),key=lambda m:m.original_id)
    with Path(path).open('w') as out:
        for m in result:
            attributes=f'gene_id "{m.gene}"; transcript_id "{m.original_id}";'
            left=min(a for a,b in m.intervals);right=max(b for a,b in m.intervals)
            out.write(f'{m.chrom}\tGAPMS\ttranscript\t{left+1}\t{right}\t.\t{m.strand}\t.\t{attributes}\n')
            for a,b in sorted(m.intervals):out.write(f'{m.chrom}\tGAPMS\texon\t{a+1}\t{b}\t.\t{m.strand}\t.\t{attributes}\n')
    write_table(provenance,['original_transcript_id','original_gene_id','transcript_id','gene_id','chrom','strand','exon_blocks'],
                sorted(records,key=lambda r:(r['transcript_id'],r['original_transcript_id'])))
    return result


def reconstruct(args, reference, lengths, commands):
    out = args.output
    directory = out/'transcripts'; directory.mkdir(exist_ok=True)
    products = out/'orfs'; products.mkdir(exist_ok=True)
    required = ['stringtie', 'TransDecoder.LongOrfs', 'TransDecoder.Predict']
    missing = [name for name in required if not shutil.which(name)]
    if missing:raise RuntimeError('BAM mode requires: '+', '.join(missing))
    transcripts = directory/'assembled_transcripts.gtf'
    marker = directory/'reconstruction.complete.json'
    raw_predictions=out/'work/rna_prediction';raw_predictions.mkdir(parents=True,exist_ok=True)
    projected = raw_predictions/'all_predicted_orfs.gtf'
    if args.resume and marker.exists() and projected.exists():
        return projected
    run_command(['stringtie', args.bam, '-p', args.threads, '-o', transcripts], directory,
                out/'logs/stringtie.log', commands)
    models = load_models(transcripts, 'RNA_transcripts', reference, lengths, args.genetic_code, allow_unstranded=True)
    shutil.copy2(transcripts,directory/'stringtie_original.gtf')
    models=canonicalize_transcripts(models,transcripts,directory/'transcript_identifiers.tsv')
    write_table(directory/'transcript_orientation.tsv', ['transcript_id','input_strand','extraction_orientation'],
                (dict(transcript_id=m.original_id,input_strand=m.strand,extraction_orientation='-' if m.strand=='-' else '+') for m in models))
    fasta = directory/'assembled_transcripts.fa'
    with fasta.open('w') as stream:
        for model in models:
            dna = join_blocks(reference, model.chrom, model.strand, model.intervals)
            stream.write(f'>{model.original_id}\n{dna}\n')
    if args.genetic_code != 1:
        raise ValueError('TransDecoder reconstruction currently supports genetic code 1 only; mapping/GTF-only modes accept other codes')
    # TransDecoder PWM code seeds its RNG, but hash traversal also affects
    # training/start refinement. Fix Perl ordering only for these child tools.
    td_env=dict(os.environ, PERL_HASH_SEED='0', PERL_PERTURB_KEYS='0')
    run_command(['TransDecoder.LongOrfs', '-t', fasta, '-m', args.min_orf_aa], directory,
                out/'logs/transdecoder_longorfs.log', commands, env=td_env)
    run_command(['TransDecoder.Predict', '-t', fasta, '--cpu', args.threads], directory,
                out/'logs/transdecoder_predict.log', commands, env=td_env)
    from .orf_projection import convert_transdecoder_to_genome_gtf
    convert_transdecoder_to_genome_gtf(str(fasta)+'.transdecoder.gff3', transcripts, projected)
    from .common import write_json
    write_json(marker, {'completed':True,'projected_gtf':str(projected)})
    return projected
