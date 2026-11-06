"""Reference/RNA-free splice-flank prior learned from molecular-qualified consensus.
Only discounts a positive protein-quality bonus; never supplies an admission vote.
"""
from pathlib import Path
from collections import Counter,defaultdict
import fcntl,hashlib,json,math,os
import numpy as np
from Bio.Seq import Seq
from .local_molecular import frame_at
from .direct_molecular import DirectMolecular
from .provenance import sha256,digest_object

BASES='ACGT';BINDEX={b:i for i,b in enumerate(BASES)};K=28

def context(genome,m,a,b):
    length=genome.lengths[m.contig]
    if m.strand=='+':
        if a<3 or a+6>length or b<20 or b+3>length:return None
        donor=genome.fetch(m.contig,a-3,a+6);acceptor=genome.fetch(m.contig,b-20,b+3)
    else:
        if b<6 or b+3>length or a<3 or a+20>length:return None
        donor=str(Seq(genome.fetch(m.contig,b-6,b+3)).reverse_complement())
        acceptor=str(Seq(genome.fetch(m.contig,a-3,a+20)).reverse_complement())
    family=donor[3:5]+acceptor[18:20]
    if family not in ['GTAG','GCAG','ATAC']:return None
    flanks=donor[:3]+donor[5:]+acceptor[:18]+acceptor[20:]
    if len(flanks)!=K or set(flanks)-set(BASES):return None
    phase=(frame_at(m,a-1 if m.strand=='+' else b)+1)%3
    key='\t'.join(map(str,(m.contig,m.strand,a,b,phase)))
    return key,family+':'+str(phase),flanks

def background(genome,genome_sha,cache):
    path=cache/('splice_background_'+genome_sha+'.json')
    with (cache/'.splice_background.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if path.exists():return json.loads(path.read_text())['probabilities']
        counts=Counter()
        for contig,length in genome.lengths.items():
            for left in range(0,length,1000000):
                seq=genome.fetch(contig,left,min(length,left+1000000))
                counts.update({b:seq.count(b) for b in BASES})
        total=sum(counts.values());p=[(counts[b]+1)/(total+4) for b in BASES]
        path.write_text(json.dumps(dict(genome_sha256=genome_sha,counts=dict(counts),probabilities=p,reference_used=False))+'\n');return p

class SplicePrior:
    def __init__(self,record):
        self.record=record;self.bg=np.array(record['background']);self.groups=record['groups'];self.training=set(record['training_keys']);self.memo={}
    def score(self,genome,m,a,b):
        data=context(genome,m,a,b)
        if data is None:return 1.,dict(observed=False,reason='unknown_or_noncanonical_context')
        key,group,flanks=data
        if key in self.memo:return self.memo[key]
        record=self.groups.get(group)
        if not record:return 1.,dict(observed=False,reason='no_qualified_training_sites',group=group)
        counts=np.array(record['counts'],dtype=np.float64);n=record['n'];x=np.array([BINDEX[b] for b in flanks]);leave_out=key in self.training
        if leave_out:counts[np.arange(K),x]-=1;n-=1
        assert (counts>=0).all()
        if n<=0:return 1.,dict(observed=False,reason='no_training_sites_after_leave_one_out',group=group)
        p=(counts+4*self.bg[None,:])/(n+4.)
        lr=float(np.log(p[np.arange(K),x]/self.bg[x]).sum())
        compatibility=1/(1+math.exp(-lr)) if lr>=0 else math.exp(lr)/(1+math.exp(lr))
        reliability=n/(n+4.*K);factor=1-reliability*(1-compatibility)
        assert 0<=factor<=1
        answer=factor,dict(observed=True,group=group,training_sites=n,leave_one_site_out=leave_out,log_likelihood_ratio=lr,compatibility_index=compatibility,reliability=reliability)
        self.memo[key]=answer;return answer

def train_or_load(originals,source_models,genome,fingerprints,cache,threshold=.437721461):
    cache=Path(cache);cache.mkdir(parents=True,exist_ok=True)
    # RNA/peptides are deliberately absent from this identity and training data.
    identity=digest_object(dict(genome=fingerprints[str(genome.path)],originals=[(c.uid,sorted(c.sources),c.q,c.complete,c.canonical) for c in sorted(originals,key=lambda c:c.uid)],code=sha256(__file__),minimum_sources=min(2,len(source_models)),threshold=threshold))
    path=cache/('intrinsic_splice_'+identity+'.json')
    with (cache/'.intrinsic_splice.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if path.exists():
            record=json.loads(path.read_text());assert record['identity']==identity;return SplicePrior(record),path
        min_sources=min(2,len(source_models));sites={};support=defaultdict(set)
        for c in originals:
            if not c.sources or c.q is None or c.q<threshold or not c.complete or not c.canonical:continue
            for a,b in c.model.introns:
                data=context(genome,c.model,a,b)
                if data is None:continue
                key,group,flanks=data;sites[key]=(group,flanks);support[key].update(c.sources)
        groups={};training=[]
        for key in sorted(sites):
            if len(support[key])<min_sources:continue
            group,flanks=sites[key];record=groups.setdefault(group,dict(n=0,counts=np.zeros((K,4),dtype=np.int64)))
            record['n']+=1;record['counts'][np.arange(K),[BINDEX[b] for b in flanks]]+=1;training.append(key)
        for record in groups.values():record['counts']=record['counts'].tolist()
        record=dict(identity=identity,reference_used=False,RNA_training=False,minimum_sources=min_sources,protein_threshold=threshold,background=background(genome,fingerprints[str(genome.path)],cache),positions=K,pseudocount_mass=4,groups=groups,training_keys=training,implementation_sha256=sha256(__file__))
        tmp=path.with_name(path.name+'.tmp.'+str(os.getpid()));tmp.write_text(json.dumps(record,separators=(',',':'))+'\n');os.replace(tmp,path)
        return SplicePrior(record),path

def attach(candidates,prior,genome):
    counts=Counter();all_factors=[]
    for c in candidates:
        factors=[]
        for a,b in c.model.introns:
            factor,detail=prior.score(genome,c.model,a,b);factors.append(factor);all_factors.append(factor);counts['measured' if detail['observed'] else 'abstained']+=1
        c.graph_evidence['intrinsic_splice_factors']=factors
    return dict(site_evaluations=dict(counts),mean_factor=float(np.mean(all_factors)) if all_factors else None,training_sites=len(prior.training),groups={k:v['n'] for k,v in prior.groups.items()},reference_used=False,RNA_training=False)
