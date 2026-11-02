"""RNA transcript reconstruction with read-only assembly access."""
from pathlib import Path
import shutil

import pysam
from gapms._peplocus.sequence import join_blocks
from .common import run_command
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


def reconstruct(args, reference, lengths, commands):
    out = args.output
    directory = out/'transcripts'; directory.mkdir(exist_ok=True)
    products = out/'orfs'; products.mkdir(exist_ok=True)
    required = ['stringtie', 'TransDecoder.LongOrfs', 'TransDecoder.Predict']
    missing = [name for name in required if not shutil.which(name)]
    if missing:raise RuntimeError('BAM mode requires: '+', '.join(missing))
    transcripts = directory/'assembled_transcripts.gtf'
    marker = directory/'reconstruction.complete.json'
    projected = products/'candidate_orfs.gtf'
    if args.resume and marker.exists() and projected.exists():
        return projected
    run_command(['stringtie', args.bam, '-p', args.threads, '-o', transcripts], directory,
                out/'logs/stringtie.log', commands)
    models = load_models(transcripts, 'RNA_transcripts', reference, lengths, args.genetic_code, allow_unstranded=True)
    from .common import write_table
    write_table(directory/'transcript_orientation.tsv', ['transcript_id','input_strand','extraction_orientation'],
                (dict(transcript_id=m.original_id,input_strand=m.strand,extraction_orientation='-' if m.strand=='-' else '+') for m in models))
    fasta = directory/'assembled_transcripts.fa'
    with fasta.open('w') as stream:
        for model in models:
            dna = join_blocks(reference, model.chrom, model.strand, model.intervals)
            stream.write(f'>{model.original_id}\n{dna}\n')
    if args.genetic_code != 1:
        raise ValueError('TransDecoder reconstruction currently supports genetic code 1 only; mapping/GTF-only modes accept other codes')
    run_command(['TransDecoder.LongOrfs', '-t', fasta, '-m', args.min_orf_aa], directory,
                out/'logs/transdecoder_longorfs.log', commands)
    run_command(['TransDecoder.Predict', '-t', fasta, '--cpu', args.threads], directory,
                out/'logs/transdecoder_predict.log', commands)
    from gapms.bam_search import convert_transdecoder_to_genome_gtf
    convert_transdecoder_to_genome_gtf(str(fasta)+'.transdecoder.gff3', transcripts, projected)
    from .common import write_json
    write_json(marker, {'completed':True,'projected_gtf':str(projected)})
    return projected
