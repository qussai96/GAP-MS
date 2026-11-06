"""Coding/noncoding competition following the documented GeneSensus convention.

Independent implementation, with phase-qualified coding/splice votes and
GAPMS's conservative peptide budget. The fixed 0.75 coefficient is inherited
from the GeneSensus/EVM convention; no reference features or fitted weights.
"""
from collections import defaultdict
from bisect import bisect_left,bisect_right
from .consensus_rna import union


class Integral:
    def __init__(self,values):
        self.blocks=union(values);self.starts=[a for a,b in self.blocks];self.ends=[b for a,b in self.blocks]
    def intersect(self,regions):
        n=0
        for a,b in regions:
            for x,y in self.blocks[bisect_right(self.ends,a):bisect_left(self.starts,b)]:n+=max(0,min(b,y)-max(a,x))
        return n


class StateContext:
    def __init__(self,candidates):
        from .consensus import offset,junctions,merge_blocks
        coding=defaultdict(list);spans=defaultdict(list);introns=defaultdict(list);oriented=defaultdict(list);exact=defaultdict(set);self.sources=set();self.exon_contexts=set()
        for c in candidates:
            if not c.valid or c.origin=='graph':continue
            for m in c.members:
                if m.issues:continue
                bb=merge_blocks(m.coding);ordered=bb if m.strand=='+' else tuple(reversed(bb));used=0
                for i,(a,b) in enumerate(ordered):
                    role='single' if len(ordered)==1 else 'first' if i==0 else 'last' if i==len(ordered)-1 else 'internal'
                    self.exon_contexts.add((m.chrom,m.strand,role,(a,b),used%3));used+=b-a
                if m.source=='RNA_ORF':continue
                self.sources.add(m.source);spans[m.source,m.chrom].append((bb[0][0],bb[-1][1]))
                for a,b in bb:
                    phase=(offset(bb,m.strand,a)+(-a if m.strand=='+' else a))%3
                    coding[m.source,m.chrom,m.strand,phase].append((a,b))
                for a,b in junctions(bb):
                    introns[m.source,m.chrom].append((a,b));oriented[m.source,m.chrom,m.strand].append((a,b))
                    phase=(offset(bb,m.strand,a-1 if m.strand=='+' else b)+1)%3
                    exact[m.chrom,m.strand,a,b,phase].add(m.source)
        self.sources=sorted(self.sources);self.coding={k:Integral(v) for k,v in coding.items()};self.spans={k:Integral(v) for k,v in spans.items()};self.introns={k:Integral(v) for k,v in introns.items()};self.oriented={k:Integral(v) for k,v in oriented.items()};self.exact=exact
    @staticmethod
    def integral(index,key,regions):
        return index[key].intersect(regions) if key in index else 0

    def features(self,c,rna):
        from .consensus import offset,junctions
        m=c.model;bb=c.blocks;length=sum(b-a for a,b in bb);span=[(bb[0][0],bb[-1][1])];extent=span[0][1]-span[0][0];jj=junctions(bb)
        coding=splice=intergenic=absent=intronic=opposite=0
        for source in self.sources:
            for a,b in bb:
                phase=(offset(bb,m.strand,a)+(-a if m.strand=='+' else a))%3
                coding+=self.integral(self.coding,(source,m.chrom,m.strand,phase),[(a,b)])
            intergenic+=extent-self.integral(self.spans,(source,m.chrom),span)
            absent+=length-self.integral(self.spans,(source,m.chrom),bb)
            intronic+=self.integral(self.introns,(source,m.chrom),bb)
            opposite+=self.integral(self.oriented,(source,m.chrom,'-' if m.strand=='+' else '+'),jj)
        for a,b in jj:
            phase=(offset(bb,m.strand,a-1 if m.strand=='+' else b)+1)%3
            splice+=(b-a)*len(self.exact.get((m.chrom,m.strand,a,b,phase),()))
        coverage=union(rna.coverage.get((m.chrom,m.strand),[])+rna.coverage.get((m.chrom,'.'),[]))
        rna_coding=Integral(coverage).intersect(bb) if coverage else 0
        rna_splice=sum((b-a)*n/(n+3.) for a,b in jj for n in [rna.count(m.chrom,m.strand,a,b)] if n>=3)
        peptide=c.features['peptide_coding_nt']+c.features['peptide_splice_nt']
        positive=coding+splice+peptide+rna_coding+rna_splice
        noncoding=absent+intronic+opposite+splice
        admitted=positive+1e-9>=.75*noncoding
        utility=positive-.75*intergenic
        # Neutral existing source hypotheses stay available without peptide evidence.
        if abs(utility)<=1e-9 and c.sources:utility=.001
        return dict(source_coding_nt=coding,source_splice_nt=splice,whole_span_intergenic_nt=intergenic,
            absent_coding_nt=absent,intronic_contradiction_nt=intronic,opposite_intronic_nt=opposite,
            admission_positive_nt=positive,equivalent_noncoding_nt=noncoding,structurally_admitted=admitted,
            rna_coding_nt=rna_coding,rna_splice_nt=rna_splice,source_utility=(coding+splice-.75*intergenic)/1000.,
            utility=max(0.,utility/1000.) if admitted else 0.)

    def supports_exons(self,c,obs):
        ordered=c.blocks if c.model.strand=='+' else tuple(reversed(c.blocks));used=0
        for i,iv in enumerate(ordered):
            role='single' if len(ordered)==1 else 'first' if i==0 else 'last' if i==len(ordered)-1 else 'internal'
            token=(c.model.chrom,c.model.strand,role,iv,used%3);used+=iv[1]-iv[0]
            if token in self.exon_contexts:continue
            if any(obs.strong(k) and iv in k[3][1:-1] for k in c.features['hits']):continue
            return False
        return True


def conservative_repair(locus,chosen,obs,rna):
    """One-for-one splice repairs preserve incumbent unique sequences and RNA."""
    from .consensus import overlap_nt,junctions,same_frame
    changes=[];chosen=list(chosen);changed=True
    while changed:
        changed=False
        for c in sorted(locus,key=lambda c:(-len(c.features['unique']),c.id)):
            if any(p.id==c.id for p in chosen) or not c.valid or c.reasons or not(c.model.complete_start and c.model.complete_stop):continue
            old=[p for p in chosen if overlap_nt(c.blocks,p.blocks)>0]
            if len(old)!=1:continue
            p=old[0]
            if not p.features['unique']<=c.features['unique'] or not same_frame(c.blocks,p.blocks,c.model.strand):continue
            diagnostic=[];split=False;newedges=set(junctions(c.blocks))-set(junctions(p.blocks))
            for key in c.features['hits']:
                if not obs.strong(key) or key in p.features['hit_set']:continue
                # Evidence must lie within the incumbent locus; no terminal-extension assertion.
                if key[3][0][0]<p.blocks[0][0] or key[3][-1][1]>p.blocks[-1][1]:continue
                diagnostic.append(key[0]);split|=bool(newedges.intersection(junctions(key[3])) and min(b-a for a,b in key[3])>=9)
            nonnested={s for s in diagnostic if not any(s!=t and s in t for t in diagnostic)}
            if not(newedges and (split or len(nonnested)>=2)):continue
            oldrna={iv for iv in junctions(p.blocks) if rna.count(p.model.chrom,p.model.strand,*iv)>=3}
            if not oldrna<=set(junctions(c.blocks)):continue
            chosen=[x for x in chosen if x.id!=p.id]+[c];changes.append(dict(candidate=c.id,incumbent=p.id,diagnostic_peptides=sorted(set(diagnostic))));changed=True;break
    return chosen,changes
