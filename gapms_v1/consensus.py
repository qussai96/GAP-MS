"""Bounded, evidence-guided consensus of supplied coding annotations.

Inspired by GeneSensus 0.6.0a3's strand-aware graph and ambiguity budget.
Scores are transparent utilities, never probabilities. Missing peptides are
unknown; original source hypotheses survive in the audit even if not selected.
"""
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from functools import cached_property
import json
import time
from urllib.parse import quote

from gapms._peplocus.sequence import canonical, join_blocks, translate
from .common import LOG, stable_id, write_json, write_table
from .mapping import locus_key
from .models import Model


def merge_blocks(blocks):
    result=[]
    for a,b in sorted(blocks):
        if result and a==result[-1][1]:result[-1]=(result[-1][0],b)
        else:result.append((a,b))
    return tuple(result)


def model_key(m):
    return m.chrom,m.strand,merge_blocks(m.coding)


def offset(blocks,strand,position):
    total=0
    for a,b in (blocks if strand=='+' else reversed(blocks)):
        if a<=position<b:return total+(position-a if strand=='+' else b-1-position)
        total+=b-a
    return None


def peptide_offset(blocks,strand,peptide_blocks):
    """Exact local path AND codon offset; overlap is insufficient."""
    pp=sorted(peptide_blocks);pos=pp[0][0] if strand=='+' else pp[-1][1]-1
    start=offset(blocks,strand,pos)
    if start is None or start%3:return None
    from gapms._peplocus.sequence import project
    ordered=blocks if strand=='+' else list(reversed(blocks))
    projected=project(start,start+sum(b-a for a,b in pp),ordered,strand)
    return start if merge_blocks(projected)==merge_blocks(pp) else None


def overlap_nt(a,b):
    i=j=n=0
    while i<len(a) and j<len(b):
        n+=max(0,min(a[i][1],b[j][1])-max(a[i][0],b[j][0]))
        if a[i][1]<b[j][1]:i+=1
        else:j+=1
    return n


def same_frame(a,b,strand):
    for x,y in a:
        for u,v in b:
            p=max(x,u)
            if p<min(y,v) and offset(a,strand,p)%3!=offset(b,strand,p)%3:return False
    return True


def canonical_junction(reference,chrom,strand,a,b):
    if b-a<20:return False
    left=reference.fetch(chrom,a,a+2).upper();right=reference.fetch(chrom,b-2,b).upper()
    if strand=='-':
        from Bio.Seq import Seq
        left,right=str(Seq(right).reverse_complement()),str(Seq(left).reverse_complement())
    return (left,right) in {('GT','AG'),('GC','AG'),('AT','AC')}


def junctions(blocks):
    return tuple((a[1],b[0]) for a,b in zip(blocks,blocks[1:]))


@dataclass
class Candidate:
    id: str
    model: Model
    members: list=field(default_factory=list)
    origin: str='source'
    sources: tuple=()
    valid: bool=True
    reasons: list=field(default_factory=list)
    features: dict=field(default_factory=dict)
    status: str='unresolved'
    parent: str=''

    @cached_property
    def blocks(self):return merge_blocks(self.model.coding)


def locus_groups(candidates):
    grouped=defaultdict(list)
    for c in candidates:grouped[c.model.chrom,c.model.strand].append(c)
    for key,values in sorted(grouped.items()):
        current=[];end=-1
        for c in sorted(values,key=lambda c:(c.blocks[0][0],c.blocks[-1][1],c.id)):
            if current and c.blocks[0][0]>=end:
                yield current;current=[];end=-1
            current.append(c);end=max(end,c.blocks[-1][1])
        if current:yield current


def make_candidate(chrom,strand,blocks,reference,code,origin='graph',members=()):
    blocks=merge_blocks(blocks);ordered=list(blocks if strand=='+' else reversed(blocks))
    uid=stable_id('C',chrom,strand,blocks)
    dna=join_blocks(reference,chrom,strand,ordered).upper();aa=translate(dna,code)
    m=Model(uid,'CONSENSUS_GRAPH',uid,uid,chrom,strand,exons=list(blocks),cds=list(blocks),coding=ordered,
            protein=aa[:-1] if aa.endswith('*') else aa,complete_start=dna.startswith('ATG'),complete_stop=len(dna)%3==0 and aa.endswith('*'))
    if '*' in m.protein:m.issues.append('internal_stop')
    if any(x not in 'ACGT' for x in dna):m.issues.append('unknown_coding_bases')
    if len(dna)%3:m.issues.append('partial_terminal_codon')
    inherited={x for source in members for x in source.issues}
    reasons=sorted(set(m.issues)|inherited)
    valid=not reasons
    if origin=='graph' and not(m.complete_start and m.complete_stop):valid=False;reasons.append('graph_ORF_incomplete')
    if origin=='graph' and not all(canonical_junction(reference,chrom,strand,a,b) for a,b in junctions(blocks)):
        valid=False;reasons.append('graph_noncanonical_splice')
    return Candidate(uid,m,list(members),origin,tuple(sorted({x.source for x in members if x.source!='RNA_ORF'})),valid,reasons)


class Observations:
    """Full mapping universe, including paths outside candidate annotations."""
    def __init__(self,store,metadata):
        self.paths={};self.bins=defaultdict(set);self.counts=store.counts();self.quality={}
        for p,item in metadata.items():self.quality.setdefault(canonical(p),set()).add(item['quality'])
        for seq,ch,st,encoded in store.db.execute('SELECT DISTINCT canonical,chrom,strand,blocks FROM hits'):
            bb=tuple(map(tuple,json.loads(encoded)));key=(seq,ch,st,bb);self.paths[key]=bb
            for a,b in bb:
                for k in range(a//1000,(b-1)//1000+1):self.bins[ch,st,k].add(key)

    def attach(self,c):
        possible=set();m=c.model
        for a,b in c.blocks:
            for k in range(a//1000,(b-1)//1000+1):possible.update(self.bins.get((m.chrom,m.strand,k),()))
        found=[]
        for key in sorted(possible):
            seq,ch,st,bb=key
            if bb[0][0]<c.blocks[0][0] or bb[-1][1]>c.blocks[-1][1]:continue
            location=peptide_offset(c.blocks,st,bb)
            if location is not None:
                # Query was genome-translated at import/projection; recheck against this model.
                start=location//3
                if canonical(m.protein[start:start+len(seq)])==seq:found.append(key)
        return found

    def strong(self,key):
        seq=key[0];q=self.quality.get(seq,set())
        return self.counts[seq]==1 and len(seq)>=9 and ('qvalue_fail' not in q or 'qvalue_pass' in q)

    def features(self,c):
        hits=self.attach(c);events=defaultdict(float);split=defaultdict(lambda:defaultdict(float));unique=set()
        for key in hits:
            seq,ch,st,bb=key;q=self.quality.get(seq,set())
            if 'qvalue_fail' in q and 'qvalue_pass' not in q:continue
            w=1/self.counts[seq]
            for a,b in bb:events[a]+=w;events[b]-=w
            jj=junctions(bb)
            for iv in jj:
                for a,b in bb:split[iv][a]+=w/len(jj);split[iv][b]-=w/len(jj)
            if self.strong(key):unique.add(seq)
        def integrate(ev):
            value=active=0.;previous=None
            for p,d in sorted(ev.items()):
                if previous is not None:value+=(p-previous)*min(1.,max(0.,active))
                active+=d;previous=p
            return value
        coverage=integrate(events);splice=sum(min(b-a,integrate(split.get((a,b),{}))) for a,b in junctions(c.blocks))
        return dict(hits=hits,hit_set=set(hits),unique=unique,peptide_coding_nt=coverage,peptide_splice_nt=splice)


def source_features(c,locus):
    coding=splice=0.
    for source in {s for d in locus for s in d.sources}:
        same=[m for d in locus for m in d.members if m.source==source and not m.issues]
        # One vote/base/source, independent of isoform count.
        ev=defaultdict(int)
        for m in same:
            for a,b in merge_blocks(m.coding):
                for u,v in c.blocks:
                    x,y=max(a,u),min(b,v)
                    if x<y and offset(c.blocks,c.model.strand,x)%3==offset(merge_blocks(m.coding),m.strand,x)%3:ev[x]+=1;ev[y]-=1
        active=0;prev=None
        for p,d in sorted(ev.items()):
            if prev is not None and active:coding+=p-prev
            active+=d;prev=p
        for a,b in junctions(c.blocks):
            upstream=a-1 if c.model.strand=='+' else b
            phase=offset(c.blocks,c.model.strand,upstream)
            if any((a,b) in junctions(merge_blocks(m.coding)) and offset(merge_blocks(m.coding),m.strand,upstream)%3==phase%3 for m in same):splice+=1
    return coding/1000.+.2*splice


def reconstruct_graph(locus,reference,observations,args):
    """Two bounded searches: structural and unique-peptide guided, union kept."""
    report=dict(locus_id=stable_id('L',sorted(c.id for c in locus)),nodes=0,paths_pruned=0,candidates=0,skipped='',phase_rejections=0)
    if len(locus)<2:return [],report
    chrom,strand=locus[0].model.chrom,locus[0].model.strand;direction=1 if strand=='+' else -1
    bounds=sorted({p for c in locus for iv in c.blocks for p in iv});report['nodes']=len(bounds)
    starts=sorted({c.blocks[0][0] if strand=='+' else c.blocks[-1][1] for c in locus},reverse=strand=='-')
    ends={c.blocks[-1][1] if strand=='+' else c.blocks[0][0] for c in locus}
    if len(bounds)>args.consensus_max_nodes or len(starts)>args.consensus_max_starts:
        report['skipped']='node_or_start_limit';return [],report
    graph=defaultdict(list)
    for a,b in zip(bounds,bounds[1:]):
        contexts=[c for c in locus if any(x<=a and b<=y for x,y in c.blocks)]
        if contexts:
            ss={s for c in contexts for s in c.sources};phases={offset(c.blocks,strand,a if strand=='+' else b-1)%3 for c in contexts}
            start,end=(a,b) if strand=='+' else (b,a)
            graph[start].append((end,'coding',a,b,phases,len(ss)))
    jjs={iv for c in locus for iv in junctions(c.blocks)}
    for a,b in sorted(jjs):
        if not canonical_junction(reference,chrom,strand,a,b):continue
        start,end=(a,b) if strand=='+' else (b,a)
        ss={s for c in locus if (a,b) in junctions(c.blocks) for s in c.sources}
        contexts=[c for c in locus if (a,b) in junctions(c.blocks)]
        phases={(offset(c.blocks,strand,a-1 if strand=='+' else b)+1)%3 for c in contexts}
        graph[start].append((end,'splice',a,b,phases,len(ss)))
    # Never splice arbitrary new exon pairs: graph edges are supplied model paths.
    known={c.blocks for c in locus};found={}
    possible=set()
    for c in locus:
        for a,b in c.blocks:
            for k in range(a//1000,(b-1)//1000+1):possible.update(observations.bins.get((chrom,strand,k),()))
    hints=[k for k in sorted(possible) if observations.strong(k)]
    for guided in [False,True]:
        if guided and not hints:continue
        states=defaultdict(list)
        for s in starts:states[s].append((0.,(),0,s,frozenset()))
        for pos in sorted(bounds,reverse=strand=='-'):
            grouped=defaultdict(dict)
            for state in states.pop(pos,[]):
                score,bb,n,start,seen=state;grouped[start,n%3][bb]=max(state,grouped[start,n%3].get(bb,state))
            current=[]
            for key,values in sorted(grouped.items()):
                ordered=sorted(values.values(),key=lambda x:(-x[0],x[1]));report['paths_pruned']+=max(0,len(ordered)-args.consensus_beam);current.extend(ordered[:args.consensus_beam])
            for score,bb,n,start,seen in current:
                terminal=(bb[-1][1] if strand=='+' else bb[0][0]) if bb else None
                if pos in ends and n%3==0 and bb and terminal==pos and bb not in known and bb not in found:
                    c=make_candidate(chrom,strand,bb,reference,args.genetic_code)
                    if c.valid:found[bb]=(score,c)
                for end,kind,a,b,phases,strength in graph.get(pos,()):
                    if kind=='splice':
                        if not bb or terminal!=pos:continue
                        if n%3 not in phases:report['phase_rejections']+=1;continue
                        newbb=bb;newn=n;bonus=.2*strength
                    else:
                        if n%3 not in phases:report['phase_rejections']+=1;continue
                        newbb=merge_blocks(bb+((a,b),));newn=n+b-a;bonus=(b-a)*strength/1000.
                    newseen=seen
                    if guided and kind=='coding':
                        matched={key[0] for key in hints if key[0] not in seen and peptide_offset(newbb,strand,key[3]) is not None}
                        if matched:
                            bonus+=sum(3*len(s)/1000.+.2*max((len(k[3])-1 for k in hints if k[0]==s),default=0) for s in matched);newseen=seen|matched
                    states[end].append((score+bonus,newbb,newn,start,newseen))
    ordered=sorted(found.values(),key=lambda x:(-x[0],x[1].id));report['candidates_before_limit']=len(ordered);report['paths_pruned']+=max(0,len(ordered)-args.consensus_max_candidates)
    chosen=[c for _,c in ordered[:args.consensus_max_candidates]];report['candidates']=len(chosen)
    return chosen,report


def prepare(args,models,store,metadata,reference):
    started=time.monotonic();groups=defaultdict(list)
    supported={row[0] for row in store.db.execute('SELECT DISTINCT model FROM associations')}
    for m in models:
        if m.source=='RNA_ORF' and m.id not in supported:continue
        if m.coding:groups[model_key(m)].append(m)
    candidates=[]
    for (ch,st,bb),members in sorted(groups.items()):
        # Conflicting source records do not poison an identical valid structure.
        preferred=sorted(members,key=lambda m:(bool(m.issues),m.source,m.id))[0]
        c=make_candidate(ch,st,bb,reference,args.genetic_code,'source',[preferred]);c.members=members;c.sources=tuple(sorted({m.source for m in members if m.source!='RNA_ORF'}))
        c.origin='rna' if not c.sources else 'source';candidates.append(c)
    obs=Observations(store,metadata);graphs=[];reports=[]
    for locus in locus_groups(candidates):
        found,report=reconstruct_graph([c for c in locus if c.valid],reference,obs,args);graphs.extend(found);reports.append(report)
    LOG.info('Consensus graph: %s source structures, %s bounded graph candidates',len(candidates),len(graphs))
    return candidates+graphs,graphs,reports,round(time.monotonic()-started,3)


def event_support(candidate,primary,obs,rna,rna_geometry=None):
    """Distinguishing local evidence; shared peptides cannot admit an isoform."""
    new=set(junctions(candidate.blocks))-set(junctions(primary.blocks))
    diagnostic=[]
    for key in candidate.features['hits']:
        if not obs.strong(key):continue
        bb=key[3]
        if key in primary.features['hit_set']:continue
        # Cross a changed junction, or code in a genuinely alternative exon region.
        if new.intersection(junctions(bb)):diagnostic.append(key[0])
    return sorted(set(diagnostic)),bool(rna_geometry and rna_geometry.supports_change(candidate,primary))


def fusion_reason(c,locus):
    """Do not select a model spanning separate source genes without connectivity proof."""
    bysource=defaultdict(list)
    for d in locus:
        if d.id==c.id:continue
        for m in d.members:
            if m.source=='RNA_ORF':continue
            if overlap_nt(c.blocks,d.blocks)>0:bysource[m.source].append((m.gene,d.blocks))
    for values in bysource.values():
        for i,(g,a) in enumerate(values):
            for h,b in values[i+1:]:
                if g!=h and overlap_nt(a,b)==0:return 'spans_separate_source_genes'
    return ''


def select_primary(locus,args):
    """Maximum utility compatible set; exact small loci, bounded large loci."""
    valid=[c for c in locus if c.valid and not c.reasons];conflicts={c.id:{d.id for d in valid if d.id!=c.id and overlap_nt(c.blocks,d.blocks)>0} for c in valid}
    order=sorted(valid,key=lambda c:(-c.features['utility'],c.id));byid={c.id:c for c in valid}
    if len(valid)<=18:
        best=(-1.,());states=0
        def visit(remaining,chosen,value):
            nonlocal best,states
            states+=1
            if value+sum(byid[x].features['utility'] for x in remaining)<best[0]-1e-9:return
            if not remaining:
                key=tuple(sorted(chosen))
                if value>best[0]+1e-9 or (abs(value-best[0])<=1e-9 and key<best[1]):best=value,key
                return
            x=remaining[0];tail=remaining[1:]
            visit([y for y in tail if y not in conflicts[x]],chosen+[x],value+byid[x].features['utility']);visit(tail,chosen,value)
        visit([c.id for c in order],[],0.)
        return [byid[x] for x in best[1]],'exact',states
    chosen=[]
    for c in order:
        if not any(d.id in conflicts[c.id] for d in chosen):chosen.append(c)
    return chosen,'greedy_bound',len(order)


def finish(args,candidates,store,metadata,reference,rna,graph_reports,graph_seconds):
    started=time.monotonic();obs=Observations(store,metadata);decisions=[];primaries=[];alternatives=[];selection_reports=[]
    from .consensus_rna import RNAEvidence
    rna_geometry=RNAEvidence.from_bam(args.bam,list(locus_groups(candidates)),args.min_mapq,args.min_read_anchor,getattr(args,'rna_library','unstranded'),args.output/'work/consensus_rna.bam.csi')
    for i,c in enumerate(candidates,1):
        if i%10000==0:LOG.info('Consensus evidence attached: %s/%s candidates',i,len(candidates))
        c.features=obs.features(c)
        c.features['rna_supported_junctions']=sum(rna_geometry.count(c.model.chrom,c.model.strand,a,b)>=3 for a,b in junctions(c.blocks))
    for locus_number,locus in enumerate(locus_groups(candidates),1):
        if locus_number%5000==0:LOG.info('Consensus selection: %s regions processed',locus_number)
        for c in locus:
            guard=fusion_reason(c,locus)
            if guard:c.reasons.append(guard)
            f=c.features;f['source_utility']=source_features(c,locus)
            # All terms are explicitly heuristic; peptide bases do not inflate to full genes.
            f['utility']=max(.001,f['source_utility']+(f['peptide_coding_nt']+f['peptide_splice_nt'])/1000.+.2*f['rna_supported_junctions'])
            if c.origin=='graph':
                # Novel combinations need local strong evidence; preserve all existing direct observations.
                hits={k[0] for k in f['hits'] if obs.strong(k)}
                split=any(obs.strong(k) and len(k[3])>1 and min(b-a for a,b in k[3])>=9 for k in f['hits'])
                nonnested=[s for s in hits if not any(s!=t and s in t for t in hits)]
                if not(split or len(nonnested)>=2 or f['rna_supported_junctions']):c.reasons.append('graph_lacks_strong_translation_or_RNA_support')
                for d in locus:
                    if d.origin=='graph' or not overlap_nt(c.blocks,d.blocks):continue
                    if not d.features['unique']<=f['unique']:c.reasons.append('graph_loses_source_unique_peptides')
                    oldrna={iv for iv in junctions(d.blocks) if rna_geometry.count(d.model.chrom,d.model.strand,*iv)>=3}
                    if not oldrna<=set(junctions(c.blocks)):c.reasons.append('graph_loses_supported_RNA_junction')
                # A graph recombination without distinguishing evidence remains provisional.
                source_competitors=[d for d in locus if d.origin!='graph' and overlap_nt(c.blocks,d.blocks)]
                if source_competitors and not any(event_support(c,d,obs,rna,rna_geometry)[0] or event_support(c,d,obs,rna,rna_geometry)[1] for d in source_competitors):c.reasons.append('graph_no_distinguishing_splice_evidence')
        chosen,solver,states=select_primary(locus,args);selection_reports.append(dict(locus_id=stable_id('L',sorted(c.id for c in locus)),solver=solver,states=states,candidates=len(locus)))
        for c in chosen:c.status='primary';c.parent=stable_id('G',c.id);primaries.append(c)
        # Retain only full source-predicted alternative splice structures with event-specific support.
        for c in sorted(locus,key=lambda c:(-c.features['utility'],c.id)):
            if any(c.id==p.id for p in chosen):continue
            overlaps=[p for p in chosen if overlap_nt(c.blocks,p.blocks)>0]
            if len(overlaps)==1 and c.origin!='graph' and c.valid and not c.reasons:
                p=overlaps[0];shared=overlap_nt(c.blocks,p.blocks)
                diagnostic,rna_ok=event_support(c,p,obs,rna,rna_geometry)
                strictshare=shared>=.5*sum(b-a for a,b in c.blocks) and shared>=.5*sum(b-a for a,b in p.blocks)
                if (junctions(c.blocks)!=junctions(p.blocks) and strictshare and same_frame(c.blocks,p.blocks,c.model.strand)
                    and (set(c.blocks)&set(p.blocks) or set(junctions(c.blocks))&set(junctions(p.blocks)))
                    and c.model.complete_start and c.model.complete_stop
                    and all(canonical_junction(reference,c.model.chrom,c.model.strand,a,b) for a,b in junctions(c.blocks))
                    and (diagnostic or rna_ok) and sum(x.parent==p.parent for x in alternatives)<2):
                    c.status='supported_alternative';c.parent=p.parent;alternatives.append(c)
                    c.features['alternative_diagnostic_peptides']=diagnostic;c.features['alternative_RNA_support']=rna_ok
            if c.reasons or not c.valid:c.status='rejected'
        for c in locus:
            competitors=[p.id for p in chosen if overlap_nt(c.blocks,p.blocks)>0 and p.id!=c.id]
            competing=[p for p in locus if p.id!=c.id and p.valid and overlap_nt(c.blocks,p.blocks)>0]
            distinguishing=[k[0] for k in c.features['hits'] if obs.strong(k) and competing
                and any(k not in p.features['hit_set'] for p in competing)]
            decisions.append(dict(candidate_id=c.id,status=c.status,selected_gene=c.parent,origin=c.origin,sources=c.sources,
                source_models=[m.id for m in c.members],chrom=c.model.chrom,strand=c.model.strand,cds_blocks=c.blocks,
                complete_start=c.model.complete_start,complete_stop=c.model.complete_stop,competitors=competitors,
                distinct_peptides=len({k[0] for k in c.features['hits']}),unique_peptides=sorted(c.features['unique']),
                distinguishing_peptides=sorted(set(distinguishing)),peptide_coding_nt=round(c.features['peptide_coding_nt'],4),
                peptide_splice_nt=round(c.features['peptide_splice_nt'],4),source_utility=round(c.features['source_utility'],6),
                rna_supported_junctions=c.features['rna_supported_junctions'],utility=round(c.features['utility'],6),
                reasons=sorted(set(c.reasons)) or ['selected_by_available_evidence_and_structural_utility' if c.status=='primary' else 'event_specific_peptide_or_RNA_support' if c.status=='supported_alternative' else 'insufficient_evidence_to_emit_additional_isoform'],
                evidence_status='distinguishing_unique_peptide' if distinguishing else 'unique_regional_support' if c.features['unique'] else 'regional_or_structural_only',
                quality='unknown' if any(obs.quality.get(s,set())=={'unknown'} for s in distinguishing) else 'see_peptide_identifications',
                boundary_status='genomic_ORF_boundaries' if c.model.complete_start and c.model.complete_stop else 'partial_source_boundaries_preserved'))
    directory=args.output/'consensus';directory.mkdir(exist_ok=True)
    selected=primaries+alternatives
    write_annotation(selected,directory/'annotation.gff3');write_annotation(primaries,directory/'primary_annotation.gff3');write_annotation(alternatives,directory/'additional_isoforms.gff3')
    write_annotation([c for c in candidates if c.status=='unresolved'],directory/'unresolved_alternatives.gff3',unresolved=True)
    with (directory/'proteins.faa').open('w') as f:
        for c in sorted(selected,key=lambda c:c.id):f.write(f'>{c.id}\n{c.model.protein}\n')
    write_table(directory/'decisions.tsv',list(decisions[0]) if decisions else ['candidate_id','status'],decisions)
    write_table(directory/'source_model_mapping.tsv',['source','original_transcript_id','source_gene','source_model_id','candidate_id','status','consensus_gene_id'],
                (dict(source=m.source,original_transcript_id=m.original_id,source_gene=m.gene,source_model_id=m.id,candidate_id=c.id,status=c.status,consensus_gene_id=c.parent) for c in candidates for m in c.members))
    write_table(directory/'graph_report.tsv',sorted({k for r in graph_reports for k in r}) or ['locus_id'],graph_reports)
    write_table(directory/'selection_report.tsv',['locus_id','solver','states','candidates'],selection_reports)
    stats=dict(policy='automatic_multi_annotation',primary_models=len(primaries),additional_isoforms=len(alternatives),selected_models=len(selected),
               candidates=len(candidates),status_counts=dict(Counter(c.status for c in candidates)),graph_candidates=sum(c.origin=='graph' for c in candidates),
               rna_qualified_reads=rna_geometry.reads,graph_selected=sum(c.origin=='graph' for c in selected),graph_seconds=graph_seconds,selection_seconds=round(time.monotonic()-started,3),
               graph_paths_pruned=sum(r['paths_pruned'] for r in graph_reports),graph_loci_skipped=sum(bool(r['skipped']) for r in graph_reports),
               bounded_selection_loci=sum(r['solver']!='exact' for r in selection_reports),reference_used=False,
               limitations=['Utilities are heuristics, not transcript correctness probabilities.','Global uniqueness is conditional on searched genome and splice-path bounds.',
                 'No peptide observations is unknown evidence, not a gene-removal rule.','Overlapping coding paths define conflicts; nested and overlapping same-strand genes can remain unresolved.',
                 'Graph recombinations require distinguishing splice evidence; remaining paths stay in diagnostics.','Partial source models are retained as partial when selected; graph models require complete ORFs.',
                 'Secondary isoforms need a distinguishing unique splice peptide or qualified RNA support for all new junctions, >=95% new-coding-base coverage and whole-retained-intron spanning observations where applicable.',
                 'Short RNA observations do not phase full transcripts; predicted starts/stops are not experimentally established.'])
    write_json(directory/'summary.json',stats);LOG.info('Consensus selected %s primaries and %s alternatives',len(primaries),len(alternatives))
    return stats


def write_annotation(candidates,path,unresolved=False):
    """Coding-only GFF3 with consistent parents, phases and genome-derived IDs."""
    genes=defaultdict(list)
    for c in candidates:genes[c.parent or stable_id('U',c.id)].append(c)
    with path.open('w') as f:
        f.write('##gff-version 3\n')
        for gid,cc in sorted(genes.items(),key=lambda x:(x[1][0].model.chrom,min(c.blocks[0][0] for c in x[1]),x[0])):
            chrom,strand=cc[0].model.chrom,cc[0].model.strand;a=min(c.blocks[0][0] for c in cc);b=max(c.blocks[-1][1] for c in cc)
            f.write(f'{chrom}\tGAPMS\tgene\t{a+1}\t{b}\t.\t{strand}\t.\tID={gid}\n')
            for c in sorted(cc,key=lambda c:c.id):
                attrs=f'ID={c.id};Parent={gid};status={c.status};origin={c.origin};sources={quote(",".join(c.sources),safe="")};partial={str(not(c.model.complete_start and c.model.complete_stop)).lower()}'
                f.write(f'{chrom}\tGAPMS\tmRNA\t{c.blocks[0][0]+1}\t{c.blocks[-1][1]}\t.\t{strand}\t.\t{attrs}\n')
                cumulative=0;phases={}
                for x,y in (c.blocks if strand=='+' else reversed(c.blocks)):phases[x,y]=(-cumulative)%3;cumulative+=y-x
                for i,(x,y) in enumerate(c.blocks,1):
                    f.write(f'{chrom}\tGAPMS\texon\t{x+1}\t{y}\t.\t{strand}\t.\tID={c.id}.exon{i};Parent={c.id}\n')
                    f.write(f'{chrom}\tGAPMS\tCDS\t{x+1}\t{y}\t.\t{strand}\t{phases[x,y]}\tParent={c.id}\n')
