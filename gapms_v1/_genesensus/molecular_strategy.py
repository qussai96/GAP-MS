"""Optional original-protein diagnostics and self-supervised genomic splice context.
No protein encoder, reference annotation, or target-trained selector is required.
"""
from pathlib import Path
import csv,gzip,hashlib,math

LOCALIZABLE={'chimeric','splice_artifact','terminal_abnormal'}
ALIASES={'spurious_orf':'cryptic_orf','spurious_orfs':'cryptic_orf','spuriousorf':'cryptic_orf','spuriousorfs':'cryptic_orf','splice':'splice_artifact','terminal':'terminal_abnormal'}


def read_diagnostics(path):
    opener=gzip.open if str(path).endswith('.gz') else open;rows={}
    with opener(path,'rt') as f:
        reader=csv.DictReader(f,delimiter='\t')
        required={'protein_id','diagnostic_class','diagnostic_class_probability'}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):raise ValueError('Diagnostics require protein_id, diagnostic_class and diagnostic_class_probability')
        for row in reader:
            ident=row['protein_id'];category=row['diagnostic_class'].strip().lower().replace(' ','_').replace('-','_')
            category=ALIASES.get(category,category)
            if not ident or category not in LOCALIZABLE|{'cryptic_orf',''}:raise ValueError('Invalid diagnostic ID/class: '+ident)
            probability=float(row['diagnostic_class_probability']) if category else None
            if category and (not math.isfinite(probability) or not 0<=probability<=1):raise ValueError('Invalid diagnostic posterior: '+ident)
            if row.get('binary_call')=='positive' and category:raise ValueError('Binary-positive proteins must not carry negative diagnostic classes')
            value=(category,probability)
            if ident in rows and rows[ident]!=value:raise ValueError('Conflicting duplicate diagnostic ID: '+ident)
            rows[ident]=value
    if not rows:raise ValueError('Empty diagnostic table')
    return rows


def prepare(candidates,source_models,config,genome,fingerprints,cache,support):
    p=config['parameters'];strategy=p['molecular_strategy']
    report=dict(strategy=strategy,enabled=False,reference_used=False,reference_trained=False,trained=False,new_protein_encoding=False,new_protein_localization=False)
    if strategy=='source_context':return candidates,source_models,None,report,[],[]
    originals=[c for c in candidates if c.sources]
    if not any(c.q is not None for c in originals):
        report['reason']='no_original_protein_scores';return candidates,source_models,None,report,[],[]
    threshold=p['localizer_gate'];assert threshold is not None
    tables={};sequence_hash={};observations={};matched=0;missing=0;source_matches={}
    for source,models in source_models.items():
        spec=config['predictions'][source];path=spec.get('diagnostics')
        if path and path not in tables:tables[path]=read_diagnostics(path)
        table=tables.get(path,{});source_matched=0
        for model in models:
            value=table.get(model.identifier)
            if value is None and table:
                if model.key not in sequence_hash:
                    aa=genome.assess(model)[2].rstrip('*');sequence_hash[model.key]=hashlib.sha256(aa.encode()).hexdigest()
                value=table.get(sequence_hash[model.key])
            observations[source,model.identifier]=value
            matched+=value is not None;missing+=value is None;source_matched+=value is not None
        source_matches[source]=dict(provided=bool(path),matched=source_matched,unobserved=len(models)-source_matched)
        if path and models and not source_matched:raise ValueError("No diagnostic identifiers or sequence hashes match source: "+source)
    labels={};allowed_tracks=set();eligible_profiles=0;abstained_profiles=0
    for c in originals:
        values=[observations.get((source,ident)) for source,ids in c.original_ids.items() for ident in ids]
        observed=[v for v in values if v is not None and v[0]]
        # A whole-model veto requires agreement of the supplied diagnostic calls.
        if observed and len(observed)==len(values) and len({v[0] for v in observed})==1:
            labels[c.uid]=(observed[0][0],min(v[1] for v in observed))
        else:labels[c.uid]=('',0.)
        retained=[]
        for track in c.risks:
            m=track.model;value=observations.get((m.source,m.identifier))
            if value and value[0] in LOCALIZABLE and track.q is not None and track.q<threshold:
                retained.append(track);allowed_tracks.add(id(track));eligible_profiles+=1
            else:abstained_profiles+=1
        c.risks=retained
    for c in candidates:
        if not c.sources:c.risks=[r for r in c.risks if id(getattr(r,'track',r)) in allowed_tracks]
    from .splice_prior import train_or_load
    from .source_splice import IntrinsicSourceMolecular
    prior,prior_path=train_or_load(originals,source_models,genome,fingerprints,cache,threshold=threshold)
    events=[];vetoed=[];orf_report=dict(vetoed_originals=0)
    if strategy=='splice_orf':
        from .orf_diagnostics import classify
        rejected,orf_report,events=classify(originals,labels,support,p,'conflict_orf_veto')
        vetoed=[c.model for c in originals if c.uid in rejected];keys={m.key for m in vetoed}
        source_models={source:[m for m in models if m.key not in keys] for source,models in source_models.items()}
        candidates=[c for c in candidates if c.uid not in rejected]
    def factory(local_originals,parameters):return IntrinsicSourceMolecular(local_originals,parameters,prior,genome)
    report.update(enabled=True,trained=bool(prior.training),training_scope='genomic splice-flank frequencies from input-predictor consensus qualified by frozen protein scores; no RNA or reference labels',training_sites=len(prior.training),training_groups={k:v['n'] for k,v in prior.groups.items()},prior_cache=str(prior_path),diagnostic_sources=source_matches,diagnostic_matches=matched,diagnostic_unobserved=missing,eligible_original_profiles=eligible_profiles,abstained_original_profiles=abstained_profiles,orf_diagnostic=orf_report)
    return candidates,source_models,factory,report,events,vetoed
