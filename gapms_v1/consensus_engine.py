"""GAPMS adapter to the frozen GeneSensus 0.6.0a3 inference engine.

The vendored algorithm is unchanged. This module normalizes GAPMS inputs,
shares its peptide mapping space and exposes the engine's original audit files.
"""
from collections import defaultdict
from dataclasses import dataclass
import csv
import gzip
import json
from pathlib import Path
import shutil
import time

from peplocus.sequence import canonical, join_blocks, translate
from .common import LOG, stable_id, write_json, write_table
from .models import Model
from ._genesensus.config import load_config, input_paths
from ._genesensus.genome import Genome
from ._genesensus.loading import load_sources
from ._genesensus.graph_cache import graph_candidates
from ._genesensus.pipeline import _rna_data, run as engine_run
from ._genesensus.peptides import read_bed
from ._genesensus.provenance import sha256
from ._genesensus.support import Support
from ._genesensus.models import Transcript


@dataclass
class GraphModel:
    model: Model


def write_evidence(store, metadata, path):
    """Keep policy-eligible placements, including off-model ones, of each eligible I/L class."""
    quality = defaultdict(set)
    for peptide, row in metadata.items():
        quality[canonical(peptide)].add(row['quality'])
    excluded = {p for p, qq in quality.items() if qq == {'qvalue_fail'}}
    rows = 0
    with path.open('w') as f:
        query = 'SELECT DISTINCT canonical,chrom,strand,blocks FROM active_hits ORDER BY chrom,strand,canonical,blocks'
        for peptide, chrom, strand, encoded in store.db.execute(query):
            if peptide in excluded:
                continue
            bb = json.loads(encoded); a, b = bb[0][0], bb[-1][1]
            fields = [chrom,a,b,peptide,0,strand,a,b,0,len(bb),
                      ','.join(str(y-x) for x,y in bb), ','.join(str(x-a) for x,y in bb)]
            f.write('\t'.join(map(str, fields))+'\n'); rows += 1
    return dict(placements=rows, excluded_failed_identification_classes=len(excluded))


def normalize_sources(args, models, directory):
    """Preserve each parsed source CDS, initial phase, exons and gene family.

    Internal phases are written consistently with the initial phase already used
    by the GAPMS mapper. Original issues/identities are kept in a separate audit.
    RNA-reconstructed ORFs are candidates, not additional predictor votes.
    """
    groups = defaultdict(list)
    for m in models:
        if m.source != 'RNA_ORF' and m.cds:
            groups[m.source].append(m)
    predictions = {}; audit = []
    for number, (source, mm) in enumerate(sorted(groups.items()), 1):
        name = f'prediction_{number:02d}'
        path = directory / (name+'.gtf')
        predictions[name] = {'annotation': str(path)}
        with path.open('w') as f:
            for m in sorted(mm, key=lambda x:x.id):
                # Safe, source-scoped IDs avoid quoting and multi-parent ambiguities.
                at = f'gene_id "{stable_id("G", source, m.chrom, m.strand, m.gene)}"; transcript_id "{m.id}";'
                ee = m.exons or m.cds
                f.write(f'{m.chrom}\tGAPMS\ttranscript\t{min(a for a,b in ee)+1}\t{max(b for a,b in ee)}\t.\t{m.strand}\t.\t{at}\n')
                for a,b in ee:
                    f.write(f'{m.chrom}\tGAPMS\texon\t{a+1}\t{b}\t.\t{m.strand}\t.\t{at}\n')
                first = m.cds[0 if m.strand == '+' else -1]
                initial = int(m.phases[first]) if m.phases[first] != '.' else 0
                used = -initial; phases = {}
                for a,b in (m.cds if m.strand == '+' else reversed(m.cds)):
                    phases[a,b] = (-used)%3; used += b-a
                for a,b in m.cds:
                    f.write(f'{m.chrom}\tGAPMS\tCDS\t{a+1}\t{b}\t.\t{m.strand}\t{phases[a,b]}\t{at}\n')
                audit.append(dict(predictor=name,source=source,source_model=m.id,original_transcript_id=m.original_id,
                    original_gene_id=m.gene,chrom=m.chrom,strand=m.strand,initial_phase=initial,
                    original_issues=m.issues,engine_candidate=Transcript(
                        m.chrom,m.strand,tuple(m.cds),name,m.id,initial,tuple(ee),m.gene).uid))
    if len(predictions) != len(args.gtf):
        raise ValueError('Each consensus annotation must contain at least one coding transcript')
    write_table(directory/'source_models.tsv',list(audit[0]),audit)
    return predictions, audit


def prepare(args, models, store, metadata, reference):
    if args.genetic_code != 1:
        raise ValueError('Multi-annotation GeneSensus consensus currently requires --genetic-code 1')
    started = time.monotonic()
    directory = args.output/'work/consensus_inputs'; directory.mkdir(parents=True,exist_ok=True)
    predictions, audit = normalize_sources(args, models, directory)
    bed = directory/'peptides.bed'; evidence = write_evidence(store,metadata,bed)
    config = dict(genome=str(Path(reference.filename.decode() if isinstance(reference.filename,bytes) else reference.filename).resolve()),
        predictions=predictions,peptides={'bed':str(bed)},parameters=dict(
            selection_policy='state_molecular_guided',molecular_strategy='source_context',
            peptide_strategy='graph_contrast',compress_diagnostics=True,
            graph_beam=args.consensus_beam,graph_candidates=args.consensus_max_candidates,
            graph_max_nodes=args.consensus_max_nodes,graph_max_starts=args.consensus_max_starts,threads=args.threads))
    if args.bam:
        # Preserve the frozen engine's RNA qualification defaults. GAPMS mapper
        # thresholds remain independent and are recorded in the run manifest.
        config['rna'] = dict(bam=str(getattr(args,'consensus_bam',None) or args.bam),library=args.rna_library)
    config_path=directory/'config.json'; write_json(config_path,config)
    config=load_config(config_path)
    cache=args.output/'work/genesensus_cache';cache.mkdir(exist_ok=True)
    fingerprints={path:sha256(path) for path in input_paths(config)}
    graphs=[]
    with Genome(config['genome']) as genome:
        sources,pool,phase,molecular,source_cache=load_sources(config,genome,fingerprints,cache)
        rna,_=_rna_data(config,cache,args.threads,genome,fingerprints)
        support=Support(rna,(),config['parameters']['minimum_rna_templates'])
        support.peptide_evidence,_=read_bed(bed,keep_ambiguous=True)
        found,reports,_=graph_candidates(pool,config,fingerprints,support,genome,cache)
        for transcript,parents in found:
            if transcript.key in pool:continue
            complete,canonical_splice,protein=genome.assess(transcript)
            coding=list(transcript.cds if transcript.strand=='+' else reversed(transcript.cds))
            trim=transcript.phase; trimmed=[]
            for a,b in coding:
                n=min(trim,b-a);trim-=n
                if transcript.strand=='+':a+=n
                else:b-=n
                if a<b:trimmed.append((a,b))
            m=Model('GS.g'+transcript.uid+'.t1','consensus_graph',transcript.identifier,transcript.gene,
                transcript.contig,transcript.strand,exons=list(transcript.exons),cds=list(transcript.cds),
                protein=protein.rstrip('*'),coding=trimmed,complete_start=complete,complete_stop=complete)
            graphs.append(GraphModel(m))
    context=dict(config_path=config_path,bed=bed,cache=cache,audit=audit,preparation_evidence=evidence,graph_ids={c.model.id for c in graphs})
    LOG.info('Frozen GeneSensus graph: %s unique new paths for shared peptide association',len(graphs))
    return context,graphs,reports,round(time.monotonic()-started,3)


def feedback_model(identifier, row, reference):
    """A graph candidate first seen after selection, with sequence-derived completeness.

    Uses the same start/stop rules as load_models; completeness is never assumed
    from the candidate's origin.
    """
    bb=[tuple(x) for x in json.loads(row['CDS_0based_halfopen'])]
    dna=join_blocks(reference,row['contig'],row['strand'],bb)
    aa=translate(dna,1)
    protein=aa[:-1] if aa.endswith('*') else aa
    issues=(['internal_stop'] if '*' in protein else [])+(['partial_terminal_codon'] if len(dna)%3 else [])
    return Model(identifier,'consensus_graph',identifier,identifier,row['contig'],row['strand'],
        exons=bb,cds=bb,protein=protein,coding=bb if row['strand']=='+' else list(reversed(bb)),issues=issues,
        complete_start=dna[:3].upper()=='ATG',complete_stop=len(dna)%3==0 and aa.endswith('*'))


def finish(args, context, store, metadata, reference, rna, graph_reports, graph_seconds, matcher=None):
    LOG.info('Selecting consensus with frozen GeneSensus 0.6.0a3 graph_contrast policy')
    directory=args.output/'consensus'
    feedback=[]
    for iteration in range(8):
        evidence=write_evidence(store,metadata,context['bed'])
        if directory.exists():
            shutil.rmtree(directory)  # this run's derived output; engine provenance is regenerated
        result=engine_run(context['config_path'],directory,threads=args.threads,cache=context['cache'])
        report=result.report
        selected=set();rows=[];additional_graphs=[]
        with gzip.open(directory/'candidates.tsv.gz','rt') as f:
            for row in csv.DictReader(f,delimiter='\t'):
                if row['selected']=='True':selected.add(row['candidate'])
                rows.append(row)
                identifier='GS.g'+row['candidate']+'.t1'
                if row['origin']!='graph' or identifier in context['graph_ids']:continue
                model=feedback_model(identifier,row,reference)
                additional_graphs.append(model);context['graph_ids'].add(identifier)
        before=store.db.execute('SELECT COUNT(*) FROM active_hits').fetchone()[0]
        if additional_graphs:
            search=args.output/'work'/f'consensus_feedback_{iteration+1}'
            if search.exists():shutil.rmtree(search)
            search.mkdir()
            queries=sorted(metadata)
            store.associate(additional_graphs,queries,search,reference,1,matcher=matcher)
        after=store.db.execute('SELECT COUNT(*) FROM active_hits').fetchone()[0]
        feedback.append(dict(pass_number=iteration+1,new_graph_proteins=len(additional_graphs),new_mapping_rows=after-before))
        if after==before:break
        LOG.info('Consensus mapping feedback: %s new placements; reselecting with updated ambiguity',after-before)
    else:
        raise RuntimeError('Consensus mapping feedback did not stabilize within eight passes; no completed annotation is claimed')
    isoforms=json.loads((directory/'isoforms.json').read_text())
    additional=set(isoforms['parent_by_candidate'])
    write_table(directory/'source_model_mapping.tsv',list(context['audit'][0])+['status'],
        (dict(row,status='supported_alternative' if row['engine_candidate'] in additional else 'primary' if row['engine_candidate'] in selected else 'not_selected') for row in context['audit']))
    # Native compressed candidates and decision JSONL remain the authoritative audit.
    write_table(directory/'decisions.tsv',['candidate_id','status','origin','sources','source_ids','cds_blocks','complete','canonical','peptides'],
        (dict(candidate_id=row['candidate'],status='supported_alternative' if row['candidate'] in additional else 'primary' if row['candidate'] in selected else 'not_selected',
          origin=row['origin'],sources=row['sources'],source_ids=row['source_ids'],cds_blocks=row['CDS_0based_halfopen'],
          complete=row['complete'],canonical=row['canonical'],peptides=row['peptides']) for row in rows))
    provisional=directory/'provisional.gff3'
    if provisional.exists():shutil.copyfile(provisional,directory/'unresolved_alternatives.gff3')
    else:(directory/'unresolved_alternatives.gff3').write_text('##gff-version 3\n')
    stats=dict(policy='automatic_multi_annotation',engine='GeneSensus',engine_version='0.6.0a3',
        engine_policy='state_molecular_guided/source_context/graph_contrast',reference_used=False,
        primary_models=report['primary_selected_count'],additional_isoforms=report['selected_count']-report['primary_selected_count'],
        selected_models=report['selected_count'],candidates=report['candidate_count'],graph_candidates=report['novel_graph_candidates'],
        graph_selected=report['selected_graph_candidates'],graph_seconds=graph_seconds,engine_seconds=report['timings']['total_seconds'],
        evidence=evidence,mapping_feedback=feedback,engine_report='report.json',
        limitations=['Peptide-only consensus is not guaranteed to outperform each input predictor.',
            'No protein molecular scores or residue profiles are imputed.',
            'The frozen secondary-isoform policy requires molecular corroboration or event-specific RNA; peptides alone do not admit secondary models.',
            'Unselected source structures remain in source_model_mapping.tsv and native candidate diagnostics, not as asserted new genes.',
            'RNA-reconstructed ORFs remain separate candidates and do not count as independent gene-predictor votes.'])
    write_json(directory/'summary.json',stats)
    return stats
