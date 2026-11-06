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
    write_table(root/'all_proteins_scores.tsv',fields,rows)


def publish(root,diag,args):
    summary=json.loads((diag/'reports/summary.json').read_text())
    for source,name in [('mapped_peptides.bed','peptides.bed'),('mapped_peptides.igv.bed','peptides.igv.bed')]:
        expose(diag/'mapping'/source,root/name)
    supported=diag/'annotation/supported_gene_models.gtf'
    if supported.exists():expose(supported,root/'supported_genes.gtf')
    if len(args.gtf)>1:
        expose(diag/'consensus/annotation.gff3',root/'consensus.gff3')
        expose(diag/'consensus/proteins.faa',root/'consensus_proteins.faa')
    if not ((root/'all_proteins_scores.tsv').exists() and not (diag/'work/source_models.gtf').exists()):
        score_table(root,diag,args)
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
    explanations={'supported_genes.gtf':'Supplied transcripts with locus-discriminating peptide support; their original structures are preserved. Ambiguous support is recorded in the score table.',
        'all_proteins_scores.tsv':'Assessment of supplied coding models, RNA ORFs when present, and selected new consensus graph proteins. NA means unavailable; classifier probabilities are not calculated in v1.',
        'peptides.bed':'Standard BED12 for downstream analysis; peptide sequences remain in protein orientation.',
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
    text+='<h2>Main results</h2><table><tr><th>File</th><th>Meaning</th></tr>'+main+'</table>'
    if (root/'all_proteins_scores.tsv').exists():text+='<p>'+link('diagnostics/score_columns.tsv','Score-table column definitions')+'</p>'
    text+='<h2>Changes to review</h2><p>'+str(len(events))+' proposal rows; these are not counts of distinct genes. Start/stop files also include explicitly labeled terminal-extension proposals. Uncertain boundaries remain unresolved.</p><ul>'+review_links+'</ul>'
    text+='<h2>Supplied annotations</h2><p>'+('These inputs were merged to produce the consensus.' if multi else 'The supplied annotation was assessed, not rewritten or automatically corrected.')+'</p><table><tr><th>Input path</th><th>SHA-256</th></tr>'+inp+'</table>'
    text+='<details><summary>Methods, limitations and detailed evidence</summary><ul>'+notes+'</ul><p>'+link('diagnostics/reports/summary.json','Run summary')+' · '+link('diagnostics/run_manifest.json','Input and software provenance')+' · '+link('diagnostics/','Detailed diagnostics')+'</p><p>Temporary work '+('is retained (--keep-work).' if args.keep_work else 'is removed after successful completion. Use --keep-work to retain it; failed runs retain their work.')+'</p></details></body></html>'
    temporary=root/'report.html.partial';temporary.write_text(text);temporary.replace(root/'report.html')
    write_table(diag/'score_columns.tsv',['column','meaning'],[dict(column='classifier_score',meaning='NA: the legacy classifier is not run by GAPMS v1; this is not a consensus probability'),dict(column='consensus_selection',meaning='Selection status of each source protein or new graph model; not an assessment of experimental truth'),dict(column='locus_discriminating_peptides',meaning='Distinct I/L-equivalent peptide sequences with one placement in the searched mapping universe'),dict(column='rna_supported_junctions',meaning='Source-model introns meeting mapping-stage RNA filters and --min-junction-support; distinct from GeneSensus consensus qualification'),dict(column='NA',meaning='Unavailable or not applicable; never an invented zero or probability')])
    write_json(diag/'output_layout.json',dict(layout_version=1,mode=mode,public_files=sorted(p.name for p in root.iterdir() if p.is_file()),review_files=sorted(p.name for p in review.glob('*')) if review.exists() else [],classifier_score_status='not_calculated',supported_genes_policy='existing_locus_discriminating_peptide_supported_transcripts'))
