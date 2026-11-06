"""Legacy rule-seeded iterative classification, independent of consensus voting."""
from collections import defaultdict
import contextlib
import io
import numpy as np
import pandas as pd
from sklearn.model_selection import RandomizedSearchCV, StratifiedKFold
from sklearn.base import clone
from xgboost import XGBClassifier
from gapms.extract import engineer_features
from gapms.model_utils import get_high_confident_proteins, get_low_confident_proteins, ITERATIVE_FEATURE_COLS
from .common import LOG, write_json
from .models import load_models, write_models


def features(models, hits, metadata):
    """Indexed peptide features; specificity is within this annotation pool."""
    proteins=defaultdict(set);genes=defaultdict(set);occurrences={}
    for m in models:
        pp=defaultdict(set)
        for h in hits.get(m.id,[]):
            if metadata[h['peptide']]['quality']=='qvalue_fail':continue
            pp[h['canonical']].add((h['aa_start'],h['aa_end']))
        occurrences[m.id]=pp
        for peptide in pp:
            proteins[peptide].add(m.id);genes[peptide].add((m.chrom,m.strand,m.gene))
    rows=[]
    for m in models:
        pp=occurrences[m.id];intervals=sorted({v for vv in pp.values() for v in vv})
        covered=0;end=0
        for a,b in intervals:covered+=max(0,b-max(a,end));end=max(end,b)
        ordered=m.coding  # Already in translation order on both strands.
        boundaries=[];used=0
        for i,(a,b) in enumerate(ordered[:-1]):
            used+=b-a
            following=ordered[i+1]
            if (b<following[0] if m.strand=='+' else following[1]<a):boundaries.append(used)
        def spans(v,j):return v[0]*3<j<v[1]*3
        supported=sum(any(all(spans(v,j) for v in vv) for vv in pp.values()) for j in boundaries)
        splice=sum(all(any(spans(v,j) for j in boundaries) for v in vv) for vv in pp.values())
        rows.append(dict(Protein=m.id,model_id=m.id,source=m.source,original_id=m.original_id,original_transcript_id=m.original_id,
            gene_id=m.gene,chrom=m.chrom,strand=m.strand,start=min(a for a,b in m.cds),end=max(b for a,b in m.cds),
            sequence_coverage=round(covered/max(1,len(m.protein)),2),protein_length=len(m.protein),
            mapped_peptides=len(pp),protein_specific_peptides=sum(len(proteins[p])==1 for p in pp),
            gene_specific_peptides=sum(len(genes[p])==1 for p in pp),splice_peptides=splice,
            internal_peptides=len(pp)-splice,splice_sites=len(boundaries),supported_splice_sites=supported,
            N_terminal_peptides=sum(all(a==0 or (a==1 and m.protein.startswith('M')) for a,b in vv) for vv in pp.values()),
            C_terminal_peptides=sum(all(b==len(m.protein) for a,b in vv) for vv in pp.values())))
    return engineer_features(pd.DataFrame(rows))


def classify(df, directory, threads):
    directory.mkdir(parents=True,exist_ok=True)
    with contextlib.redirect_stdout(io.StringIO()):
        high=set(get_high_confident_proteins(df).index);low=set(get_low_confident_proteins(df).index)-high
    initial=pd.Series('unlabeled',index=df.index)
    initial.loc[list(high)]='high';initial.loc[list(low)]='low'
    labels=pd.Series(np.nan,index=df.index);labels.loc[list(high)]=1;labels.loc[list(low)]=0
    df=df.copy();df['seed_label']=initial;df['classifier_score']='NA'
    audit=dict(high_seeds=len(high),low_seeds=len(low),initially_unlabeled=int((initial=='unlabeled').sum()),iterations=[],status='insufficient_seeds')
    if min(len(high),len(low))>=2:
        cv=StratifiedKFold(n_splits=min(5,len(high),len(low)),shuffle=True,random_state=42)
        # Parallelize folds, each estimator uses one CPU: respect the job allocation.
        search=RandomizedSearchCV(XGBClassifier(eval_metric='logloss',random_state=42,n_jobs=1),
            dict(n_estimators=[100,200],max_depth=[3,5,7],learning_rate=[.01,.1],subsample=[.8,1.],
                 colsample_bytree=[.8,1.],scale_pos_weight=[1,2,5]),n_iter=20,scoring='roc_auc',cv=cv,
            random_state=42,n_jobs=threads,pre_dispatch=threads)
        ix=labels.notna();search.fit(df.loc[ix,ITERATIVE_FEATURE_COLS],labels[ix])
        model=search.best_estimator_;audit.update(status='trained',best_parameters=search.best_params_,seed_cv_auroc=float(search.best_score_))
        for iteration in range(1,6):
            ix=labels.isna()
            if not ix.any():break
            probabilities=pd.Series(model.predict_proba(df.loc[ix,ITERATIVE_FEATURE_COLS])[:,1],index=df.index[ix])
            pos=probabilities.index[probabilities>=.9];neg=probabilities.index[probabilities<=.1]
            audit['iterations'].append(dict(iteration=iteration,new_positive=len(pos),new_negative=len(neg)))
            if not len(pos)+len(neg):break
            labels.loc[pos]=1;labels.loc[neg]=0;ix=labels.notna()
            model=clone(model);model.fit(df.loc[ix,ITERATIVE_FEATURE_COLS],labels[ix])
        df['classifier_score']=model.predict_proba(df[ITERATIVE_FEATURE_COLS])[:,1]
        model.save_model(directory/'model.json')
    df['iterative_label']=labels.fillna(-1).astype(int)
    positive=df['classifier_score'].map(lambda p:False if p=='NA' else float(p)>=.5)
    df['supported']=(initial=='high')|((initial=='unlabeled')&positive)
    df['confidence']=np.where(initial=='high','high',np.where(df['supported'],'moderate',np.where(initial=='low','low','unresolved' if audit['status']!='trained' else 'not_supported')))
    audit['supported']=int(df.supported.sum());write_json(directory/'summary.json',audit)
    df.to_csv(directory/'scores.tsv',sep='\t',index=False)
    LOG.info('Classifier %s: %s high seeds, %s low seeds, %s supported',audit['status'],len(high),len(low),audit['supported'])
    return df,audit


def run(args, models, store, metadata, peptides, matcher, reference, lengths):
    directory=args.output/'classification';directory.mkdir(exist_ok=True)
    pools=defaultdict(list)
    for m in models:
        if m.protein:pools[m.source].append(m)
    if len(args.gtf)>1:
        annotation=args.output/'consensus/annotation.gff3'
        has_cds=any('\tCDS\t' in line for line in annotation.open())
        selected=load_models(annotation,'CONSENSUS',reference,lengths,args.genetic_code) if has_cds else []
        for m in selected:m.id=m.original_id
        selected=[m for m in selected if m.protein]
        if selected:
            work=args.output/'work/classifier_consensus';work.mkdir(exist_ok=True)
            before=store.counts();store.associate(selected,peptides,work,reference,args.genetic_code,matcher=matcher)
            if before!=store.counts():raise RuntimeError('Classification changed the peptide mapping universe')
            pools['CONSENSUS']=selected
    hits=store.model_hits();tables=[];summaries={};supported=[]
    for source,pool in sorted(pools.items()):
        df,audit=classify(features(pool,hits,metadata),directory/source,args.threads)
        df['classification_pool']=source;tables.append(df);summaries[source]=audit
        if source=='CONSENSUS' or (len(args.gtf)<2 and source!='RNA_ORF'):
            keep=set(df.loc[df.supported,'model_id']);supported.extend(m for m in pool if m.id in keep)
        elif not args.gtf and source=='RNA_ORF':
            keep=set(df.loc[df.supported,'model_id']);supported.extend(m for m in pool if m.id in keep)
    if tables:pd.concat(tables,ignore_index=True).to_csv(directory/'scores.tsv',sep='\t',index=False)
    write_models(supported,directory/'supported_genes.gtf');write_json(directory/'summary.json',summaries)
    return summaries
