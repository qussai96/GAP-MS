"""Coding/splice versus noncoding utility, integrated with local path optimization.

Uses sparse interval integrals, capped one vote/source. No annotation reference.
EVM's published-source default coding/noncoding ratio 0.75 is an explicit fixed
compatibility convention, not a parameter fitted against this benchmark.
"""
from collections import defaultdict
from math import log1p
from bisect import bisect_left,bisect_right
from .intervals import merge,intersection_length
from .molecular import local_risk

class StateEvidence:
    def __init__(self,source_models):
        self.source_exons=set()
        coding=defaultdict(list);spans=defaultdict(list);introns=defaultdict(list);oriented=defaultdict(list);exact=defaultdict(set)
        for source,models in source_models.items():
            for m in models:
                self.source_exons.update((m.contig,m.strand,a,b) for a,b in m.cds)
                coding[source,m.contig,m.strand].extend(m.cds)
                spans[source,m.contig].append((m.cds[0][0],m.cds[-1][1]))
                introns[source,m.contig].extend(m.introns)
                oriented[source,m.contig,m.strand].extend(m.introns)
                for iv in m.introns:exact[m.contig,m.strand,iv].add(source)
        self.sources=sorted(source_models);self.exact=exact
        self.oriented_introns=self.index(oriented)
        self.coding=self.index(coding);self.spans=self.index(spans);self.introns=self.index(introns)
    @staticmethod
    def index(data):
        return {k:(iv,tuple(a for a,b in iv),tuple(b for a,b in iv)) for k,blocks in data.items() for iv in [merge(blocks)]}
    @staticmethod
    def integral(index,key,regions):
        record=index.get(key)
        if not record:return 0
        blocks,starts,ends=record;answer=0
        for a,b in regions:
            lo=bisect_right(ends,a);hi=bisect_left(starts,b)
            answer+=intersection_length(((a,b),),blocks[lo:hi])
        return answer
    def features(self,c,peptides,support,p,molecular_states=None):
        m=c.model;span=((m.cds[0][0],m.cds[-1][1]),);genomic=span[0][1]-span[0][0]
        coding=sum(self.integral(self.coding,(s,m.contig,m.strand),m.cds) for s in self.sources)
        splice=sum((b-a)*len(self.exact.get((m.contig,m.strand,(a,b)),())) for a,b in m.introns)
        intergenic=sum(genomic-self.integral(self.spans,(s,m.contig),span) for s in self.sources)
        intronic=sum(self.integral(self.introns,(s,m.contig),m.cds) for s in self.sources)
        # Contradiction over intronic sequence includes only nonmatching splice states.
        # Exact same-strand introns offset their own noncoding contribution.
        rival_intron=0
        for source in self.sources:
            covered=self.integral(self.introns,(source,m.contig),m.introns)
            shared=sum(b-a for a,b in m.introns if source in self.exact.get((m.contig,m.strand,(a,b)),()))
            rival_intron+=max(0,covered-shared)
        # Intergenic absence is relevant to candidate coding, not its declared introns.
        absent_coding=sum(m.length-self.integral(self.spans,(s,m.contig),m.cds) for s in self.sources)
        noncoding=absent_coding+intronic+rival_intron
        risk,known=local_risk(c,m.cds,p['localizer_gate']);q=c.q*(1-p['risk_strength']*(risk if known else 0)) if c.q is not None else 0.
        pep=len(c.peptides)/(len(c.peptides)+p['peptide_saturation']);rna=support.rna_score(m);protein=support.protein_score(m)
        # Bounded per-model optional channels in coding-base equivalents.
        weak_coding=weak_splice=0.
        rna_geometry=bool(support.rna_blocks)
        rna_coding=support.rna_coding_overlap(m) if rna_geometry else 0
        rna_splice=sum((b-a)*n/(n+3.) for a,b in m.introns for n in [support.rna_count(m,(a,b))] if n>=p['minimum_rna_templates']) if rna_geometry else 0.
        rna_utility=rna_coding+rna_splice if rna_geometry else m.length*rna
        if p['selection_policy'] in ('state_rna_counts','state_rna_log','state_gap','state_gap_soft','state_admitted','state_coherent','state_soft_admitted','state_splice_validated','state_molecular_guided'):
            transform = log1p if p['selection_policy']=='state_rna_log' else float
            rna_coding=sum((b-a)*transform(support.rna_exon_counts.get((m.contig,m.strand,a,b),0)) for a,b in m.cds)
            rna_splice=sum((b-a)*transform(support.rna_chain_counts.get((m.contig,m.strand,a,b),0)) for a,b in m.introns if support.rna_count(m,(a,b))>=p['minimum_rna_templates'])
            if p['selection_policy']=='state_molecular_guided':
                old_rna_coding,old_rna_splice=rna_coding,rna_splice
                # Weak positive RNA can corroborate an independently predicted
                # exact exon/junction; discovery and physical splice guards stay strict.
                rna_coding=sum((b-a)*(support.rna_exon_counts_all.get((m.contig,m.strand,a,b),0) if (m.contig,m.strand,a,b) in self.source_exons else support.rna_exon_counts.get((m.contig,m.strand,a,b),0)) for a,b in m.cds)
                rna_splice=sum((b-a)*support.rna_chain_counts.get((m.contig,m.strand,a,b),0) for a,b in m.introns if self.exact.get((m.contig,m.strand,(a,b))) or support.rna_count(m,(a,b))>=p['minimum_rna_templates'])
            if p['selection_policy']=='state_molecular_guided':weak_coding,weak_splice=rna_coding-old_rna_coding,rna_splice-old_rna_splice
            rna_utility=rna_coding+rna_splice
        # A peptide supports its translated blocks, never the full protein length.
        # Union prevents overlapping peptides and repeated sequences inflating coverage.
        peptide_blocks=merge([iv for seq in c.peptides for iv in peptides[seq].blocks])
        peptide_nt=intersection_length(m.cds,peptide_blocks)
        peptide_junctions={ (left[1],right[0]) for seq in c.peptides for left,right in zip(peptides[seq].blocks,peptides[seq].blocks[1:]) }
        peptide_splice=sum(b-a for a,b in m.introns if (a,b) in peptide_junctions)
        if p.get('peptide_strategy','legacy') != 'legacy' and hasattr(peptides,'features'):
            peptide_nt,peptide_splice,peptide_sequences=peptides.features(c,
                ambiguous=p['peptide_strategy'] != 'bounded_unique')
            pep=peptide_sequences/(peptide_sequences+p['peptide_saturation'])
        molecular_coding,molecular_splice=molecular_states.features(m) if molecular_states is not None else (m.length*q,0.)
        extra=p['molecular_weight']*(molecular_coding+molecular_splice)+m.length*p['protein_weight']*protein+p['peptide_weight']*(peptide_nt+peptide_splice)+p['rna_weight']*rna_utility
        # Two predeclared hypotheses: signed state utility vs EVM-like ratio admission.
        positive=coding+splice+extra
        opposite='-' if m.strand=='+' else '+'
        opposite_intronic=sum(self.integral(self.oriented_introns,(source,m.contig,opposite),m.introns) for source in self.sources)
        equivalent_noncoding=absent_coding+intronic+opposite_intronic+splice
        # Protein plausibility qualifies supplied states, but adds no independent
        # structural observation capable of satisfying structural admission.
        admission_positive=coding+splice+m.length*p['protein_weight']*protein+p['peptide_weight']*(peptide_nt+peptide_splice)+p['rna_weight']*rna_utility
        admitted=admission_positive>=.75*equivalent_noncoding and m.length>=p.get('minimum_coding_nt',150)
        if p['selection_policy'] in ('state_gap','state_gap_soft','state_admitted','state_coherent','state_soft_admitted','state_splice_validated','state_molecular_guided'):
            utility=positive-(.75 if p['selection_policy'] in ('state_gap_soft','state_soft_admitted','state_splice_validated','state_molecular_guided') else 1.)*intergenic
            # Preserve supplied hypotheses in exact evidence-neutral ties.
            if abs(utility)<1e-9 and c.sources:utility=1e-6
        elif p['selection_policy'] in ('state_signed','state_rna_counts','state_rna_log'):utility=positive-.75*noncoding
        else:utility=positive if positive>=.75*noncoding else 0.
        if p['selection_policy'] in ('state_admitted','state_coherent','state_soft_admitted','state_splice_validated','state_molecular_guided') and not admitted:utility=0.
        return dict(admission_positive_evidence=admission_positive,molecular_counts_for_admission=False,molecular_measured_candidate_q=c.q,molecular_source_context_coding=molecular_coding,molecular_source_context_splice=molecular_splice,rna_weak_corroboration_coding=weak_coding,rna_weak_corroboration_splice=weak_splice,equivalent_noncoding=equivalent_noncoding,coding_noncoding_ratio=positive/max(equivalent_noncoding,1e-12),admitted=admitted,coding_support=coding,exact_splice_support=splice,absent_coding=absent_coding,intronic_contradiction=intronic,rival_intron=rival_intron,noncoding=noncoding,whole_span_intergenic_opportunity_cost=intergenic,optional_base_equivalents=extra,rna_observed_coding=rna_coding,rna_observed_splice=rna_splice,molecular=(molecular_coding+molecular_splice)/max(m.length,1),peptide=pep,peptide_supported_nt=peptide_nt,peptide_supported_splice=peptide_splice,rna=rna,protein=protein,utility=max(0.,utility),unclipped_utility=utility)
