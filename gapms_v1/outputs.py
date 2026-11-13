"""Small public result set and an offline HTML entry point.

This is presentation only: native mapping, assessment and consensus decisions
remain unchanged under diagnostics. Public scientific files share disk blocks
with their native copies when hard links are available.
"""
import csv
import gzip
import hashlib
import html
import json
import os
from pathlib import Path
import shutil
from urllib.parse import quote

from .common import read_table,write_table,write_json
from .assessment import EVENT_FIELDS,MODEL_FIELDS
from .models import attributes


def records(path):
    return list(read_table(path)) if path.is_file() else []


def expose(source,destination):
    if destination.exists():destination.unlink()
    try:os.link(source,destination)
    except OSError:shutil.copyfile(source,destination)


def decoded(row,key):
    value=row.get(key,'')
    return json.loads(value) if value else []


def review_features(rows,path):
    """Emit supported local geometry; a peptide is never expanded into a CDS."""
    with path.open('w') as f:
        f.write('##gff-version 3\n')
        for e in rows:
            eid=e['event_id'];kind=e['event_type'];chrom=e['chrom'];strand=e['strand']
            at=';'.join(k+'='+quote(str(v),safe='') for k,v in dict(event_id=eid,event_type=kind,
                confidence=e['confidence'],boundary_status=e.get('boundary_status',''),
                whole_model_status=e.get('whole_model_status',''),proposal_requires_review='true',
                reference_models=e.get('reference_model_ids',''),supporting_peptides=e.get('peptides','')).items())
            if kind=='novel_orf':
                bb=decoded(e,'candidate_cds_blocks')
                if bb:
                    gene=eid+'.gene';tid=eid+'.transcript'
                    f.write(f'{chrom}\tGAPMS\tgene\t{bb[0][0]+1}\t{bb[-1][1]}\t.\t{strand}\t.\tID={gene};{at}\n')
                    f.write(f'{chrom}\tGAPMS\tmRNA\t{bb[0][0]+1}\t{bb[-1][1]}\t.\t{strand}\t.\tID={tid};Parent={gene};{at}\n')
                    phases={};used=0
                    for a,b in (bb if strand=='+' else reversed(bb)):phases[a,b]=(-used)%3;used+=b-a
                    for i,(a,b) in enumerate(bb,1):
                        f.write(f'{chrom}\tGAPMS\texon\t{a+1}\t{b}\t.\t{strand}\t.\tID={tid}.exon{i};Parent={tid}\n')
                        f.write(f'{chrom}\tGAPMS\tCDS\t{a+1}\t{b}\t.\t{strand}\t{phases[a,b]}\tParent={tid}\n')
                    continue
            features=decoded(e,'proposed_feature_blocks')
            precise=bool(features) and kind in ('alternative_start','alternative_stop') and sum(b-a for a,b in features)==3
            feature='start_codon' if precise and kind=='alternative_start' else 'stop_codon' if precise else 'peptide_supported_region'
            for i,(a,b) in enumerate(features or decoded(e,'local_peptide_blocks'),1):
                f.write(f'{chrom}\tGAPMS\t{feature}\t{a+1}\t{b}\t.\t{strand}\t.\tID={eid}.region{i};{at}\n')
            for i,(a,b) in enumerate(decoded(e,'diagnostic_junctions'),1):
                f.write(f'{chrom}\tGAPMS\tsplice_junction\t{a+1}\t{b}\t.\t{strand}\t.\tID={eid}.junction{i};{at}\n')


def rna_support(diag,minimum):
    path=diag/'transcripts/rnaseq_junctions.tsv'
    if not path.exists():return {}
    counts={}
    for row in read_table(path):
        counts[row['chrom'],int(row['start']),int(row['end']),row['strand']]=int(row['support'])
    by_model={}
    path=diag/'work/source_models.gtf'
    if not path.exists():return {}
    with path.open() as f:
        for line in f:
            if line.startswith('#'):continue
            row=line.rstrip().split('\t')
            if len(row)!=9 or row[2]!='CDS':continue
            mid=attributes(row[8])['transcript_id']
            by_model.setdefault((mid,row[0],row[6]),set()).add((int(row[3])-1,int(row[4])))
    result={}
    for (mid,ch,st),bb in by_model.items():
        bb=sorted(bb);values=[max(counts.get((ch,a[1],b[0],st),0),counts.get((ch,a[1],b[0],'.'),0)) for a,b in zip(bb,bb[1:]) if a[1]<b[0]]
        result[mid]=dict(rna_supported_junctions=sum(n>=minimum for n in values),rna_minimum_junction_support=min(values) if values else 'NA')
    return result


def score_table(root,diag,args):
    rows=records(diag/'annotation/gene_model_evidence.tsv')+records(diag/'orfs/all_predicted_orf_evidence.tsv')
    rows=[row for row in rows if row['category']!='no_CDS_annotation']
    if not rows:return
    selection={r['source_model']:r['status'] for r in records(diag/'consensus/source_model_mapping.tsv')}
    rna=rna_support(diag,args.min_junction_support) if args.bam else {}
    fields=MODEL_FIELDS+['consensus_selection','classifier_score','rna_supported_junctions','rna_minimum_junction_support']
    for row in rows:
        row.update(consensus_selection=selection.get(row['model_id'],'not_applicable' if len(args.gtf)<2 else 'not_a_predictor_vote'),classifier_score='NA',**rna.get(row['model_id'],dict(rna_supported_junctions='NA',rna_minimum_junction_support='NA')))
    # New graph proteins have no original-source table row. Native peptide count
    # counts unique sequences, not an invented total including ambiguous support.
    path=diag/'consensus/candidates.tsv.gz'
    if path.exists():
        with gzip.open(path,'rt') as f:
            for c in csv.DictReader(f,delimiter='\t'):
                if c['origin']!='graph' or c['selected']!='True':continue
                bb=json.loads(c['CDS_0based_halfopen']);row={k:'NA' for k in fields}
                row.update(model_id='GS.g'+c['candidate']+'.t1',source='consensus_graph',chrom=c['contig'],strand=c['strand'],start=bb[0][0],end=bb[-1][1],category='selected_consensus_graph',locus_discriminating_peptides=c['peptides'],consensus_selection='primary',classifier_score='NA',reasons='See native GeneSensus evidence; unknown fields are not imputed')
                rows.append(row)
    classification=records(diag/'classification/scores.tsv')
    by_id={r['model_id']:r for r in rows}
    extra=[]
    for c in classification:
        row=by_id.get(c['model_id'])
        if row is None:
            row={k:'NA' for k in fields};row.update(model_id=c['model_id'],source=c['source'],category='consensus_model',consensus_selection='primary');rows.append(row)
        for key,value in c.items():
            if key not in fields and key not in extra:extra.append(key)
            row[key]=value
    fields+=extra
    for row in rows:
        for key in extra:row.setdefault(key,'NA')
    write_table(diag/'classification/model_details.tsv',fields,rows)
    legacy_scores(root,diag,args)



def publish(root,diag,args, *, refresh=False):
    summary=json.loads((diag/'reports/summary.json').read_text())
    for source,name in [('mapped_peptides.bed','peptides.bed'),('mapped_peptides.igv.bed','peptides.igv.bed'),('tryptic_peptides.bed','peptides.tryptic.bed')]:
        expose(diag/'mapping'/source,root/name)
    supported=diag/'classification/supported_genes.gtf'
    if supported.exists():expose(supported,root/'supported_genes.gtf')
    if len(args.gtf)>1:
        expose(diag/'consensus/annotation.gff3',root/'consensus.gff3')
        expose(diag/'consensus/proteins.faa',root/'consensus_proteins.faa')
    if not refresh:score_table(root,diag,args)
    rescue={r['peptide']:r for r in records(diag/'mapping/rescue/rescue_peptides.tsv')}
    unmapped=records(diag/'mapping/unmapped_peptides.tsv')
    for row in unmapped:
        n=int(rescue.get(row['peptide'],{}).get('rescue_loci',0));row['rescue_loci']=n
        row['status']='rescue_only' if n else 'unmapped_in_searched_space'
    for row in records(diag/'mapping/excluded_peptides.tsv'):
        unmapped.append(dict(peptide_id='NA',peptide=row['peptide'],length=row['length'],status='not_searched_below_min_length',rescue_loci=0))
    write_table(root/'unmapped_peptides.tsv',['peptide_id','peptide','length','status','rescue_loci'],unmapped)
    events=records(diag/'refinement/annotation_change_proposals.tsv')
    events+=records(diag/'orfs/novel_orf_candidates.tsv')
    # The annotation-change table can already contain RNA-derived novel ORFs.
    # Export/count each event once, retaining distinct biological hypotheses.
    events=list({row['event_id']:row for row in events}.values())
    review=root/'review'
    if review.exists():shutil.rmtree(review)
    if events:
        review.mkdir();write_table(review/'candidates.tsv',EVENT_FIELDS,events)
        for name,kinds in [('alternative_start',{'alternative_start','n_terminal_extension'}),('alternative_stop',{'alternative_stop','c_terminal_extension'}),('alternative_splice',{'alternative_splice'}),('novel_orfs',{'novel_orf'})]:
            selected=[e for e in events if e['event_type'] in kinds]
            if selected:review_features(selected,review/(name+'.gff3'))
    staged=diag/'work/consensus_inputs'
    if staged.exists():
        kept=diag/'consensus/inputs';kept.mkdir(exist_ok=True)
        relocated={}
        for source in sorted(staged.iterdir()):
            if not source.is_file():continue
            destination=kept/source.name;expose(source,destination)
            relocated[str(source)]=str(destination.relative_to(diag))
        write_json(diag/'consensus/retained_input_paths.json',dict(relocated=relocated,original_assembly=str(args.assembly.resolve()),note='Native configuration and provenance retain the original execution paths; this map locates their retained annotation and peptide inputs after temporary-work cleanup.'))
    inputs=[]
    for p in args.gtf:
        p=p.resolve();h=hashlib.sha256()
        with p.open('rb') as f:
            for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
        inputs.append(dict(path=str(p),sha256=h.hexdigest()))
    write_table(diag/'input_annotations.tsv',['path','sha256'],inputs)
    multi=len(args.gtf)>1;mode='Multi-annotation consensus' if multi else 'Single-annotation assessment' if args.gtf else 'RNA and peptide mapping' if args.bam else 'Peptide mapping'
    explanations={'supported_genes.gtf':'High-confidence rule seeds plus classifier-positive initially unlabeled transcripts; consensus models in multi-annotation runs. Structures are preserved.',
        'all_proteins_scores.tsv':'Legacy-format protein features and final confidence calls for the assessed annotation. Additional classifier and source details are in diagnostics/classification/.',
        'peptides.bed':'Standard BED12 for downstream analysis; peptide sequences remain in protein orientation.',
        'peptides.tryptic.bed':'Tryptic-compatible exact placements; multiplicity and BED scores are counted independently within this subset.',
        'peptides.igv.bed':'IGV display track; negative-strand labels are reversed for visualization.',
        'unmapped_peptides.tsv':'Queries without an exact mapping, including separately labeled rescue-only and below-minimum-length sequences.',
        'consensus.gff3':'New consensus from all supplied predictions. Review proposals have not been automatically applied.',
        'consensus_proteins.faa':'Proteins corresponding to the consensus annotation.','report.log':'Progress, warnings, errors and completion status.'}
    def link(path,label=None):return '<a href="'+quote(path,safe='/')+'">'+html.escape(label or path)+'</a>'
    main=''.join('<tr><td>'+link(n)+'</td><td>'+html.escape(d)+'</td></tr>' for n,d in explanations.items() if (root/n).exists())
    review_links=''.join('<li>'+link(str(p.relative_to(root)))+'</li>' for p in sorted(review.glob('*'))) if review.exists() else '<li>No supported proposals were reported.</li>'
    inp=''.join('<tr><td><a href="'+html.escape(Path(r['path']).as_uri(),quote=True)+'">'+html.escape(r['path'])+'</a></td><td><code>'+r['sha256']+'</code></td></tr>' for r in inputs)
    m=summary['mapping'];cards=f"<p><strong>{m['mapped_peptides']:,} / {m['peptides']:,}</strong> searched peptide sequences mapped; <strong>{m['unique_peptides']:,}</strong> have one placement and <strong>{m['multiple_peptides']:,}</strong> have multiple placements.</p>"
    notes=''.join('<li>'+html.escape(x)+'</li>' for x in summary['scientific_limitations'])
    text='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>GAPMS results</title><style>body{font:16px system-ui,sans-serif;max-width:1080px;margin:40px auto;padding:0 24px;color:#24313b;line-height:1.6}h1,h2{color:#136b5b}table{border-collapse:collapse;width:100%}td,th{padding:10px;text-align:left;border-bottom:1px solid #d8e1e7;overflow-wrap:anywhere}a{color:#0867a5}code{font-size:12px}aside{padding:12px 18px;background:#eff7f4;border-left:4px solid #198870}</style><body>'
    text+='<h1>GAPMS '+html.escape(summary['version'])+'</h1><h2>'+mode+'</h2>'+cards
    text+='<aside>All supplied isoforms are checked before proposing changes. Peptide compatibility supports observed regions; it does not validate every exon, start, stop or complete transcript. Review GFF3s are proposals, not automatic corrections.</aside>'
    details=json.loads((diag/'mapping/mapping_summary.json').read_text())
    if details['engine'].get('search_scope',{}).get('splice_path_cap_reached') or summary.get('annotation_guided_search',{}).get('splice_contexts_capped'):
        text+='<aside>Splice search limits were reached. Placement uniqueness is conditional on the searched paths; inspect the diagnostic summary before interpreting proposals.</aside>'
    from .report_figure import overview
    figure=overview(root,diag,args,summary,events)
    text+=figure
    text+='<h2>Main results</h2><table><tr><th>File</th><th>Meaning</th></tr>'+main+'</table>'
    if (root/'all_proteins_scores.tsv').exists():text+='<p>'+link('diagnostics/score_columns.tsv','Score-table column definitions')+'</p>'
    specific_events=sum(e['confidence']!='ambiguous' for e in events)
    text+='<h2>Changes to review</h2><p>'+str(len(events))+' proposal rows ('+str(specific_events)+' with specific support; '+str(len(events)-specific_events)+' with ambiguous support); these are not counts of distinct genes. Start/stop files also include explicitly labeled terminal-extension proposals. Uncertain boundaries remain unresolved.</p><ul>'+review_links+'</ul>'
    text+='<h2>Supplied annotations</h2><p>'+('These inputs were merged to produce the consensus.' if multi else 'The supplied annotation was assessed, not rewritten or automatically corrected.')+'</p><table><tr><th>Input path</th><th>SHA-256</th></tr>'+inp+'</table>'
    text+='<details><summary>Methods, limitations and detailed evidence</summary><ul>'+notes+'</ul><p>'+link('diagnostics/reports/summary.json','Run summary')+' · '+link('diagnostics/run_manifest.json','Input and software provenance')+' · '+link('diagnostics/','Detailed diagnostics')+'</p><p>Temporary work '+('is retained (--keep-work).' if args.keep_work else 'is removed after successful completion. Use --keep-work to retain it; failed runs retain their work.')+'</p></details></body></html>'
    temporary=root/'report.html.partial';temporary.write_text(text);temporary.replace(root/'report.html')
    write_table(diag/'score_columns.tsv',['column','meaning'],[dict(column='classifier_score',meaning='Rule-seeded iterative XGBoost score; NA when fewer than two high or low seeds. Not a calibrated biological or consensus probability'),dict(column='consensus_selection',meaning='Selection status of each source protein or new graph model; not an assessment of experimental truth'),dict(column='locus_discriminating_peptides',meaning='Distinct I/L-equivalent peptide sequences with one placement in the searched mapping universe'),dict(column='rna_supported_junctions',meaning='Source-model introns meeting mapping-stage RNA filters and --min-junction-support; distinct from GeneSensus consensus qualification'),dict(column='NA',meaning='Unavailable or not applicable; never an invented zero or probability')])
    write_json(diag/'output_layout.json',dict(layout_version=1,mode=mode,public_files=sorted(p.name for p in root.iterdir() if p.is_file()),review_files=sorted(p.name for p in review.glob('*')) if review.exists() else [],classifier_score_status='rule_seeded_iterative',supported_genes_policy='high_seeds_plus_initially_unlabeled_classifier_positive'))


LEGACY_SCORE_FIELDS = ['Protein','supported','confidence','protein_length','isoforms','splice_sites',
    'transcript_score','sequence_coverage','mapped_peptides','N_terminal_peptides','C_terminal_peptides',
    'protein_specific_peptides','gene_specific_peptides','splice_peptides','internal_peptides','external_score',
    'mapped_peptides_per_aa','protein_specific_rate','gene_specific_rate','splice_rate','internal_rate',
    'terminal_rate','isoform_specificity','junction_evidence','has_N_term','has_C_term','has_both_term',
    'has_any_protein_specific','has_any_splice']


def assessed_scores(diag,args):
    rows=records(diag/'classification/scores.tsv')
    if len(args.gtf)>1:return [r for r in rows if r['classification_pool']=='CONSENSUS']
    if args.gtf:return [r for r in rows if r['classification_pool']!='RNA_ORF']
    return rows


def legacy_scores(root,diag,args):
    rows=assessed_scores(diag,args)
    genes={};scores={}
    for row in rows:
        key=(row['classification_pool'],row['chrom'],row['strand'],row['gene_id'])
        genes.setdefault(key,set()).add(row['Protein'])
    # Prediction scores are metadata, never invented classifier probabilities.
    for i,path in enumerate(sorted(args.gtf,key=lambda p:str(p.resolve())),1):
        source=f'A{i:02d}_{path.name}'
        opener=gzip.open if str(path).endswith('.gz') else open
        with opener(path,'rt') as f:
            for line in f:
                if line.startswith('##FASTA'):break
                if line.startswith('#'):continue
                parts=line.rstrip().split('\t')
                if len(parts)!=9 or parts[2].lower() not in {'transcript','mrna'}:continue
                at=attributes(parts[8]);tid=at.get('id') or at.get('transcript_id')
                try:score=float(parts[5])
                except ValueError:score=0.
                scores[source,parts[0],parts[6],tid]=score
    public=[]
    for row in rows:
        row=dict(row);source=row['classification_pool']
        row['isoforms']=len(genes[source,row['chrom'],row['strand'],row['gene_id']])
        row['transcript_score']=scores.get((source,row['chrom'],row['strand'],row['original_transcript_id']),0.)
        row['external_score']=0
        row['supported']='+' if row['supported']=='True' else '-'
        row['confidence']='high' if row['seed_label']=='high' else 'moderate' if row['supported']=='+' else 'no-support'
        public.append(row)
    public.sort(key=lambda r:(r['supported']!='+',r['confidence'],-float(r['mapped_peptides']),r['Protein']))
    if rows or args.gtf:write_table(root/'all_proteins_scores.tsv',LEGACY_SCORE_FIELDS,public)


def refresh_report(output):
    """Rebuild presentation from completed diagnostics without rerunning science.

    Core implementation provenance and cached score/RNA evidence are preserved.
    Publication code identity is recorded separately after a successful refresh.
    """
    from types import SimpleNamespace
    from datetime import datetime, timezone
    from . import __version__
    root=Path(output).resolve();diag=root/'diagnostics';path=diag/'run_manifest.json'
    manifest=json.loads(path.read_text())
    if manifest.get('status')!='completed':raise ValueError('report refresh requires a completed run')
    params=dict(manifest['parameters']);params['output']=root
    for key in ['assembly','peptides','bam']:
        if params.get(key) is not None:params[key]=Path(params[key])
    params['gtf']=[Path(p) for p in params.get('gtf',[])]
    params['keep_work']=bool(manifest.get('work_retained',False))
    publish(root,diag,SimpleNamespace(**params),refresh=True)
    manifest['publication']=dict(version=__version__,refreshed_utc=datetime.now(timezone.utc).isoformat(),
        implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()+Path(__file__).with_name('report_figure.py').read_bytes()).hexdigest(),
        preserved_core_implementation_sha256=manifest['implementation_sha256'],
        note='Report/event deduplication refresh only; mapping, model assessment, classifier scores and RNA evidence were preserved.')
    write_json(path,manifest)
    return manifest['publication']
