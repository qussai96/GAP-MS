"""Unified four-mode orchestration; legacy workflow code is left intact."""
from dataclasses import replace
from datetime import datetime, timezone
import hashlib
import json
import logging
import os
from pathlib import Path
import shutil
import subprocess
import time

import pysam
from gapms._peplocus.io import prepare_reference
from gapms._peplocus.junctions import Junction, extract_bam, write_junctions
from gapms._peplocus.annotation import read_annotation
from gapms._peplocus.cli import main as peplocus_main
from . import __version__
from .common import LOG, input_identity, stable_id, write_json, write_table
from .peptides import load_peptides
from .models import load_models, write_models
from .mapping import EvidenceStore
from .assessment import assess
from .refinement import refine_annotations
from . import consensus_engine as consensus
from gapms.protein_mapper import ProteinMatcher
from .rnaseq import reconstruct, validate_bam


def fingerprint():
    digest=hashlib.sha256()
    import gapms
    for root in [Path(__file__).parent,Path(gapms.__file__).parent]:
        for path in sorted(root.rglob('*.py')):
            digest.update(str(path.relative_to(root)).encode());digest.update(path.read_bytes())
    return digest.hexdigest()


def tool_versions(bam):
    result={'pysam':pysam.__version__}
    if bam:
        for name,flags in [('stringtie',['--version']),('TransDecoder.LongOrfs',['--version']),('TransDecoder.Predict',['--version'])]:
            executable=shutil.which(name)
            if not executable:raise ValueError(f'missing BAM reconstruction tool: {name}')
            try:
                process=subprocess.run([executable]+flags,capture_output=True,text=True,timeout=15)
                result[name]={'path':executable,'version_output':(process.stdout+process.stderr).strip()[:500]}
            except subprocess.TimeoutExpired:
                result[name]={'path':executable,'version_output':'version_probe_timed_out'}
    return result


def _run(args):
    started=time.monotonic()
    args.output=args.output.resolve()
    if len(args.gtf)>1 and args.genetic_code!=1:
        raise ValueError('Multi-annotation GeneSensus consensus currently requires --genetic-code 1')
    if args.bam and args.genetic_code!=1:
        raise ValueError('BAM ORF prediction currently requires --genetic-code 1')
    inputs=[args.assembly,args.peptides]+args.gtf+([args.bam] if args.bam else [])
    for path in inputs:
        if not path.is_file():raise ValueError(f'input does not exist: {path}')
    # Resolve all paths before external commands change their working directory.
    args.assembly=args.assembly.resolve();args.peptides=args.peptides.resolve()
    args.gtf=[p.resolve() for p in args.gtf];args.bam=args.bam.resolve() if args.bam else None
    if len(set(args.gtf))!=len(args.gtf):raise ValueError('same annotation supplied more than once')
    params={k:str(v) if isinstance(v,Path) else [str(x) for x in v] if k=='gtf' else v for k,v in vars(args).items() if k not in {'resume','output'}}
    identities=[input_identity(p) for p in inputs]
    run_key=stable_id('R',params,identities,fingerprint())
    manifest_path=args.output/'run_manifest.json'
    previous=None
    if args.output.exists() and any(args.output.iterdir()):
        if not args.resume or not manifest_path.exists():raise ValueError('output is nonempty; choose a new directory or use --resume for the same run')
        previous=json.loads(manifest_path.read_text())
        if previous['run_key']!=run_key:raise ValueError('resume refused: inputs, parameters or implementation differ')
        if previous['status']=='completed':
            LOG.info('Run already completed: %s',args.output);return
    args.output.mkdir(parents=True,exist_ok=True)
    for directory in ['logs','mapping','reports','work']:(args.output/directory).mkdir(exist_ok=True)
    handler=logging.FileHandler(args.output/'logs/workflow.log')
    handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
    LOG.setLevel(logging.INFO);LOG.addHandler(handler)
    if not any(isinstance(h,logging.StreamHandler) and not isinstance(h,logging.FileHandler) for h in LOG.handlers):
        stream=logging.StreamHandler();stream.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'));LOG.addHandler(stream)
    manifest=dict(version=__version__,run_key=run_key,status='running',started_utc=datetime.now(timezone.utc).isoformat(),
                  parameters=params,inputs=identities,implementation_sha256=fingerprint(),commands=previous.get('commands',[]) if previous else [],
                  host=os.uname().nodename,slurm_job_id=os.environ.get('SLURM_JOB_ID'),mode='integrated' if args.gtf and args.bam else 'annotation' if args.gtf else 'rnaseq' if args.bam else 'mapping')
    write_json(manifest_path,manifest)
    store=None
    try:
        manifest['external_tools']=tool_versions(args.bam)
        manifest['automatic_consensus']=len(args.gtf)>1
        LOG.info('GAPMS v%s %s workflow; %s threads',__version__,manifest['mode'],args.threads)
        peptides,metadata,input_stats,normalized=load_peptides(args)
        genome,index,lengths=prepare_reference(args.assembly,args.output/'work/reference')
        manifest['assembly_contigs']=len(lengths);manifest['assembly_bases']=sum(lengths.values())
        with pysam.FastaFile(genome,filepath_index=index) as reference:
            if args.bam:validate_bam(args.bam,lengths)
            models=[];source_rows=[]
            for i,path in enumerate(sorted(args.gtf,key=lambda p:str(p)),1):
                source=f'A{i:02d}_{path.name}'
                models.extend(load_models(path,source,reference,lengths,args.genetic_code))
                source_rows.append(dict(source=source,path=str(path),kind='input_annotation'))
            rna_path=None
            if args.bam:
                LOG.info('Reconstructing RNA-supported transcripts and candidate ORFs')
                rna_path=reconstruct(args,reference,lengths,manifest['commands'])
                models.extend(load_models(rna_path,'RNA_ORF',reference,lengths,args.genetic_code))
                source_rows.append(dict(source='RNA_ORF',path=str(rna_path),kind='StringTie_TransDecoder'))
            write_table(args.output/'reports/annotation_sources.tsv',['source','path','kind'],source_rows)
            rna_junctions=[];rna_stats={}
            rna_file=args.output/'transcripts/rnaseq_junctions.tsv'
            if args.bam:
                from gapms._peplocus.junctions import read_junctions
                cache=args.output/'work/rna_junction_stats.json'
                if args.resume and rna_file.exists() and cache.exists():
                    rna_junctions=read_junctions(rna_file,reference,lengths,min_support=args.min_junction_support)
                    rna_stats=json.loads(cache.read_text())
                else:
                    LOG.info('Extracting RNA splice-junction observations')
                    rna_junctions,rna_stats=extract_bam(args.bam,reference,lengths,args.min_mapq,args.min_read_anchor,args.min_junction_support)
                    write_junctions(rna_file,rna_junctions);write_json(cache,rna_stats)
            junctions={(j.chrom,j.start,j.end,j.strand):j for j in rna_junctions}
            if models:
                combined=args.output/'work/source_models.gtf';write_models(models,combined)
                annotation_junctions,annotation_stats=read_annotation(combined,reference,lengths)
                for j in annotation_junctions:
                    key=(j.chrom,j.start,j.end,j.strand)
                    old=junctions.get(key)
                    junctions[key]=replace(j,source='annotation,bam',support=old.support) if old else j
                manifest['annotation_junction_statistics']=annotation_stats
            junction_file=args.output/'work/search_junctions.tsv';write_junctions(junction_file,junctions.values())
            engine=args.output/'work/peplocus'
            if not(args.resume and (engine/'summary.json').exists()):
                if engine.exists():shutil.rmtree(engine)
                argv=['map','--genome',str(args.assembly),'--peptides',str(normalized),'--out',str(engine),'--threads',str(args.threads),
                      '--genetic-code',str(args.genetic_code),'--min-length',str(args.min_length),'--igv-display-bed',
                      '--max-junctions',str(args.max_junctions),'--max-paths-per-junction',str(args.max_paths_per_junction)]
                if junctions:argv+=['--junctions',str(junction_file)]
                if args.rescue_one_edit:argv+=['--rescue-one-edit','--rescue-min-length',str(args.rescue_min_length)]
                manifest['commands'].append({'component':'bundled_peplocus','argv':argv})
                write_json(manifest_path,manifest)
                LOG.info('Mapping peptides against the assembly and unified splice paths')
                result=peplocus_main(argv)
                if result:raise RuntimeError(f'PepLocus failed with exit {result}')
            engine_summary=json.loads((engine/'summary.json').read_text())
            args.splice_search_capped=engine_summary.get('search_scope',{}).get('splice_path_cap_reached',False)
            store_path=args.output/'work/evidence.sqlite'
            for suffix in ['', '-wal','-shm']:
                Path(str(store_path)+suffix).unlink(missing_ok=True)
            store=EvidenceStore(store_path)
            store.import_genomic(engine/'mappings.tsv')
            protein_matcher=ProteinMatcher(peptides) if models else None
            if models:
                LOG.info('Associating peptides with %s source-aware coding models',len(models))
                protein_dir=args.output/'work/protein_search'
                if protein_dir.exists():shutil.rmtree(protein_dir)
                protein_dir.mkdir()
                store.associate(models,peptides,protein_dir,reference,args.genetic_code,matcher=protein_matcher)
            consensus_stats={}
            if len(args.gtf)>1:
                LOG.info('Multiple annotation inputs: automatic consensus is enabled')
                consensus_candidates,graph_models,graph_reports,graph_seconds=consensus.prepare(args,models,store,metadata,reference)
                if graph_models:
                    graph_dir=args.output/'work/consensus_protein_search'
                    if graph_dir.exists():shutil.rmtree(graph_dir)
                    graph_dir.mkdir()
                    store.associate([c.model for c in graph_models],peptides,graph_dir,reference,args.genetic_code,matcher=protein_matcher)
            guided_events=[];refinement_stats={}
            if args.gtf:
                guided_events,refinement_stats=refine_annotations(args,models,store,metadata,reference,{(j.chrom,j.start,j.end,j.strand):j.support for j in rna_junctions})
                if refinement_stats.get('splice_contexts_capped'):args.splice_search_capped=True
            if len(args.gtf)>1:
                consensus_stats=consensus.finish(args,consensus_candidates,store,metadata,reference,
                    {(j.chrom,j.start,j.end,j.strand):j.support for j in rna_junctions},graph_reports,graph_seconds,matcher=protein_matcher)
            mapping_stats=store.export(args.output/'mapping',peptides,metadata)
            model_stats=assess(args,models,store,metadata,{(j.chrom,j.start,j.end,j.strand):j.support for j in rna_junctions},guided_events) if models else {}
        limits=json.loads((engine/'summary.json').read_text())
        write_json(args.output/'mapping/mapping_summary.json',dict(mapping=mapping_stats,engine=limits,
                   annotation_guided_search=refinement_stats,interpretation='Exact genomic placements, verified local annotation-guided paths and full-CDS projections; rescue is separate and excluded from model/event evidence.'))
        if args.rescue_one_edit:
            rescue=args.output/'mapping/rescue';rescue.mkdir(exist_ok=True)
            for path in engine.glob('rescue_*'):
                if path.is_file():shutil.copyfile(path,rescue/path.name)
        summary=dict(version=__version__,mode=manifest['mode'],consensus=consensus_stats,input=input_stats,mapping=mapping_stats,assessment=model_stats,annotation_guided_search=refinement_stats,rna_junctions=rna_stats,
                     seconds=round(time.monotonic()-started,3),scientific_limitations=[
                       'Mapping validates sequence compatibility, not peptide identification accuracy.',
                       'I/L-equivalent sequences are one evidence unit; reported multiplicity is within the searched genome and supported splice paths.',
                       'No peptide support does not imply an incorrect annotation.',
                       'RNA-supported local splice paths do not prove full-length transcript connectivity.',
                       'Novel and boundary events are candidates for annotation review, not experimentally confirmed changes.',
                       'Partial ORFs and uncertain peptide quality do not receive high-confidence candidate status.',
                       'Full unannotated ORF reconstruction requires BAM; supplied annotations also support bounded local refinement without BAM.',
                       'The legacy XGBoost classifier remains available through gapms and is not applied to novel ORFs.'])
        write_json(args.output/'reports/summary.json',summary)
        after=[input_identity(p) for p in inputs]
        if identities!=after:raise RuntimeError('an input identity changed during the run')
        manifest.update(status='core_completed',completed_utc=datetime.now(timezone.utc).isoformat(),seconds=summary['seconds'])
        write_json(manifest_path,manifest)
        LOG.info('Completed: %s/%s peptide sequences mapped; results: %s',mapping_stats['mapped_peptides'],len(peptides),args.output)
    except BaseException as error:
        manifest.update(status='failed',error=str(error),failed_utc=datetime.now(timezone.utc).isoformat())
        write_json(manifest_path,manifest)
        LOG.exception('Workflow failed')
        raise
    finally:
        if store is not None:store.db.close()
        LOG.removeHandler(handler);handler.close()


def run(args):
    """Keep the public directory small; execute native stages inside diagnostics."""
    import copy
    from .outputs import publish
    root=args.output.resolve();diag=root/'diagnostics'
    if root.exists() and any(root.iterdir()):
        if not args.resume:
            raise ValueError('output is nonempty; choose a new directory or use --resume for the same run')
        if not (diag/'run_manifest.json').exists() and any(p.name not in {'report.log','diagnostics'} for p in root.iterdir()):
            raise ValueError('resume requires an existing GAPMS run manifest')
    root.mkdir(parents=True,exist_ok=True)
    handler=logging.FileHandler(root/'report.log');handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(name)s: %(message)s'))
    logger=logging.getLogger();level=logger.level;logger.setLevel(logging.INFO);logger.addHandler(handler)
    internal=copy.copy(args);internal.output=diag
    prior=diag/'run_manifest.json'
    already_completed=prior.exists() and json.loads(prior.read_text()).get('status')=='completed'
    try:
        _run(internal)
        if args.resume and already_completed and (root/'report.html').exists():
            LOG.info('Results available: %s',root/'report.html');return
        LOG.info('Publishing the result overview and review files')
        publish(root,diag,args)
        work=diag/'work'
        if not args.keep_work and work.exists():
            try:shutil.rmtree(work)
            except OSError as error:LOG.warning('Temporary work could not be fully removed: %s',error)
        manifest_path=diag/'run_manifest.json';manifest=json.loads(manifest_path.read_text())
        manifest.update(status='completed',public_output=str(root),output_layout_version=1,work_retained=work.exists());write_json(manifest_path,manifest)
        LOG.info('Results ready: %s',root/'report.html')
    except BaseException as error:
        path=diag/'run_manifest.json'
        if path.exists():
            manifest=json.loads(path.read_text());manifest.update(status='failed',error=str(error));write_json(path,manifest)
        LOG.exception('Run failed; work and diagnostics retained')
        raise
    finally:
        logger.removeHandler(handler);handler.close();logger.setLevel(level)
