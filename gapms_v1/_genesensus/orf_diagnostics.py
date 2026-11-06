"""Conservative whole-ORF diagnostic hypotheses; reference-free, source-counted."""
from collections import defaultdict,Counter
from .graph import loci
from .intervals import overlap
from .local_molecular import frame_at


def classify(originals,diagnostics,support,p,arm):
    assert arm in ['isolated_orf_veto','conflict_orf_veto']
    bykey={c.model.key:c for c in originals};components=[]
    for models in loci([c.model for c in originals]):
        cc=[bykey[m.key] for m in models];parent=list(range(len(cc)))
        def root(i):
            while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
            return i
        for i,c in enumerate(cc):
            for j in range(i):
                if overlap(c.model.cds,cc[j].model.cds):parent[root(i)]=root(j)
        groups=defaultdict(list)
        for i,c in enumerate(cc):groups[root(i)].append(c)
        components.extend(groups.values())
    rejected=set();events=[];stats=Counter()
    for cc in components:
        sources=set().union(*(c.sources for c in cc));eligible={};protected={}
        for c in cc:
            category,posterior=diagnostics.get(c.uid,('',0.))
            eligible[c.uid]=bool(c.gene_scores and all(q<p['localizer_gate'] for q in c.gene_scores) and category=='cryptic_orf' and posterior>.5)
            reasons=[]
            if len(c.sources)>=2:reasons.append('multiple_predictors_exact_model')
            m=c.model
            for other in cc:
                if not(other.sources-c.sources):continue
                n=other.model
                if m.introns and m.introns==n.introns:reasons.append('multiple_predictors_intron_chain');break
                if not m.introns and not n.introns:
                    a,b=m.cds[0];x,y=n.cds[0]
                    if ((a<=x and y<=b) or (x<=a and b<=y)) and frame_at(m,max(a,x))==frame_at(n,max(a,x)):
                        reasons.append('multiple_predictors_nested_same_frame_coding');break
            peptide_evidence=getattr(support,'peptide_evidence',None)
            peptide_positive=(peptide_evidence.strong_translation(c)
                if peptide_evidence is not None and hasattr(peptide_evidence,'strong_translation')
                else bool(c.peptides))
            if peptide_positive:reasons.append('positive_peptide')
            if any(support.rna_count(m,iv)>=p['minimum_rna_templates'] for iv in m.introns):reasons.append('positive_exact_rna_junction')
            if any(support.rna_exon_counts_all.get((m.contig,m.strand,a,b),0)>=p['minimum_rna_templates'] for a,b in m.cds):reasons.append('positive_exact_rna_exon')
            if any((m.contig,m.strand,a,b) in getattr(support,'protein_introns',set()) for a,b in m.introns):reasons.append('positive_exact_protein_alignment_junction')
            if not m.introns and hasattr(support,'protein_score') and support.protein_score(m)==1.:reasons.append('complete_protein_alignment_coverage')
            protected[c.uid]=reasons
        all_bad=bool(len(sources)>=2 and len(cc)>=2 and all(eligible[c.uid] and not protected[c.uid] and len(c.sources)==1 for c in cc))
        for c in cc:
            veto=eligible[c.uid] and not protected[c.uid] and (len(sources)==1 or arm=='conflict_orf_veto' and all_bad)
            if veto:rejected.add(c.uid)
            stats['eligible_originals']+=eligible[c.uid];stats['protected_eligible_originals']+=bool(eligible[c.uid] and protected[c.uid]);stats['vetoed_originals']+=bool(veto)
            if eligible[c.uid]:events.append(dict(candidate=c.uid,q=c.q,diagnostic_class=diagnostics[c.uid][0],diagnostic_posterior=diagnostics[c.uid][1],sources=sorted(c.sources),component_sources=sorted(sources),component_size=len(cc),protection=protected[c.uid],veto=bool(veto),reason='isolated_single_predictor_spurious_orf' if veto and len(sources)==1 else 'all_spurious_conflicting_predictors' if veto else 'protected_or_abstained'))
    return rejected,dict(stats),events
