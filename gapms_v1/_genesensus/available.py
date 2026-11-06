"""Corroborate local predictor agreement; retain provisional genes explicitly."""
from .exchange import LocalVotes,compare_exchange
from .selection import compare,_peptide_pair,_rna_pair
from .molecular import local_risk
from .intervals import disputed,intersection_length
from functools import partial


def compare_available(c,old,peptides,support,chains,ivotes,n_sources,p,votes):
    regions=disputed(c.model,old.model,p['boundary_window_nt'])
    sn,so,known,certificate=votes.structural(c.model,old.model)
    pn,po,_,_=_peptide_pair(c,old,peptides,regions,p['peptide_saturation'])
    rn,ro=_rna_pair(c,old,support)
    cr=tuple((max(a,x),min(b,y)) for a,b in c.model.cds for x,y in regions if max(a,x)<min(b,y))
    oc=tuple((max(a,x),min(b,y)) for a,b in old.model.cds for x,y in regions if max(a,x)<min(b,y))
    an,ao=support.protein_score(c.model,cr),support.protein_score(old.model,oc)
    nr,nk=local_risk(c,regions,p['localizer_gate']);orr,ok=local_risk(old,regions,p['localizer_gate'])
    gate=p.get('confidence_gate')
    molecular_consistent=c.q is not None and old.q is not None
    corroboration=bool(pn>po or rn>ro or an>=ao+.1 or molecular_consistent or(nk and ok and orr-nr>=p['risk_advantage']))
    comparable=(c.q is not None and old.q is not None) or(nk and ok) or bool(pn or po) or abs(rn-ro)>1e-12 or abs(an-ao)>1e-12
    if known and abs(sn-so)>1e-12 and comparable:
        e=compare_exchange(c,old,peptides,support,chains,ivotes,n_sources,p,votes)
        e['decision_mode']='availability_aware_joint_exchange'
    else:
        e=compare(c,old,peptides,support,chains,ivotes,n_sources,p)
        e['decision_mode']='compatibility_without_corroborated_majority'
    e['external_channel_comparable']=bool(comparable)
    e['positive_distinguishing_evidence']=corroboration
    return e


def select_available(candidates,source_models,peptides,support,p):
    from .selection import select
    votes=LocalVotes(source_models)
    selected,events,backbone=select(candidates,source_models,peptides,support,p,
        _comparison=partial(compare_available,votes=votes))
    if p['selection_policy']!='available_tiered':return selected,events,backbone
    gate=p.get('confidence_gate')
    if gate is None:raise ValueError('Tiered selection requires an independently specified confidence_gate')
    primary=[]
    for c in selected:
        m=c.model
        # Independent predictor contradiction is unavailable with fewer than 3 sources.
        fractions=[]
        for source in votes.sources:
            blocks,starts,ends=votes.coverage.get((source,m.contig,m.strand),((),(),()))
            from bisect import bisect_left,bisect_right
            sliced=blocks[bisect_right(ends,m.cds[0][0]):bisect_left(starts,m.cds[-1][1])]
            fractions.append(intersection_length(m.cds,sliced)/m.length)
        coding=sum(fractions)/len(fractions)
        rna=any(support.rna_count(m,iv)>=p['minimum_rna_templates'] for iv in m.introns)
        protein=support.protein_score(m)>0
        consensus_provisional=(len(source_models)>=3 and coding<.5 and not(c.peptides or rna or protein))
        if consensus_provisional:
            events.append(dict(event='consensus_provisional',candidate=c.uid,incumbents=[],accepted=True,reasons=['predictor_majority_exclusion_without_positive_empirical_support'],coding_agreement=coding,interpretation='optional stricter structural-confidence view; not a false-gene label'))
        provisional=(consensus_provisional and not c.graph_path and c.q is not None and c.q<gate)
        if provisional:
            events.append(dict(event='provisional',candidate=c.uid,incumbents=[],accepted=True,
                reasons=['classifier_negative_and_predictor_majority_exclusion_without_positive_empirical_support'],
                gene_q=c.q,confidence_gate=gate,coding_agreement=coding,
                interpretation='uncertain hypothesis retained in provisional.gff3; not a biological false-gene label'))
        else:primary.append(c)
    return primary,events,backbone
