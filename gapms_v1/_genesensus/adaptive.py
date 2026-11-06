"""Experimental reference-free local evidence policy; not calibrated probabilities.

Missing peptide/RNA observations never independently delete an original model.
Peptides must distinguish the actual coding dispute; graph paths have no invented q.
"""
from .intervals import CodingIndex, disputed, intersection_length, merge
from .molecular import local_risk
from .selection import predictor_support, choose_backbone, chain_key, _structure


def peptide_dispute(c, old, peptides, regions):
    """Sequence-deduplicated, unique-placement peptides already frame-validated."""
    def side(a, b):
        seqs = sorted(s for s in a.peptides - b.peptides
                      if intersection_length(peptides[s].blocks, regions))
        blocks = merge(iv for s in seqs for iv in peptides[s].blocks)
        coding = merge((max(x, a0), min(y, b0)) for x,y in a.model.cds
                       for a0,b0 in regions if max(x,a0) < min(y,b0))
        size = sum(y-x for x,y in coding)
        coverage = intersection_length(blocks, coding)/size if size else 0.
        junctions = set(a.model.introns)-set(b.model.introns)
        bridged = sorted(iv for iv in junctions if any(
            iv in tuple((x[1], y[0]) for x,y in zip(peptides[s].blocks, peptides[s].blocks[1:]))
            for s in seqs))
        independent = len(blocks)
        strong = bool(bridged) or (len(seqs)>=2 and independent>=2 and coverage>=0.5)
        # A local short peptide can supply a small score, but cannot alone pass
        # the change gate unless it spans the differing junction.
        value = max(coverage*len(seqs)/(len(seqs)+2.),
                    len(bridged)/len(junctions) if junctions else 0.)
        return dict(sequences=seqs, coverage=coverage, bridged_junctions=bridged,
                    independent_blocks=independent, strong=strong, value=value)
    return side(c, old), side(old, c)


def compare_adaptive(c, old, peptides, support, chains, votes, n_sources, p):
    m,o=c.model,old.model
    regions=disputed(m,o,p['boundary_window_nt'])
    pn,po=peptide_dispute(c,old,peptides,regions)
    nr,nk=local_risk(c,regions,p['localizer_gate'])
    rr,rk=local_risk(old,regions,p['localizer_gate'])
    channels={}
    if c.q is not None and old.q is not None:
        channels['gene_quality']=(c.q,old.q,p['molecular_weight']*.5)
    if nk and rk:
        channels['residue_quality']=(1-nr,1-rr,p['molecular_weight']*.5)
    if pn['sequences'] or po['sequences']:
        channels['peptide']=(pn['value'],po['value'],p['peptide_weight'])
    ni=set(m.introns)-set(o.introns);oi=set(o.introns)-set(m.introns)
    rn,ro=support.rna_score(m,ni),support.rna_score(o,oi)
    if rn or ro:
        channels['rna']=(rn,ro,p['rna_weight'])
    an,ao=support.protein_score(m,regions),support.protein_score(o,regions)
    if an or ao:
        channels['protein']=(an,ao,p['protein_weight'])
    sn,so=_structure(c,n_sources,chains),_structure(old,n_sources,chains)
    if sn!=so: channels['structure']=(sn,so,p['structure_weight'])
    # Half-weights give q and residue profiles one joint molecular budget.
    # If only one is measured, it receives the full molecular budget.
    measured=[k for k in ('gene_quality','residue_quality') if k in channels]
    if len(measured)==1:
        k=measured[0];a,b,w=channels[k];channels[k]=(a,b,2*w)
    den=sum(w for a,b,w in channels.values())
    delta=sum(w*(a-b) for a,b,w in channels.values())/den if den else 0.
    molecular=(c.q is not None and old.q is not None and c.q-old.q>=p['q_advantage'])
    residue=nk and rk and rr-nr>=p['risk_advantage']
    rna=bool(ni) and rn>ro and all(support.rna_count(m,iv)>=support.min_rna_templates for iv in ni)
    protein=an>=ao+.1
    structural=(c.complete and not old.complete and
                len(chains.get(chain_key(m),()))>=2)
    distinguishing=any((molecular,residue,rna,protein,pn['strong'],structural))
    reasons=[]
    if not(c.complete and c.canonical): reasons.append('incomplete_or_noncanonical')
    if m.length<p['minimum_retained_fraction']*o.length: reasons.append('truncation')
    if not distinguishing: reasons.append('insufficient_distinguishing_evidence')
    # Never erase a unique peptide on the incumbent to satisfy another channel.
    if po['sequences']: reasons.append('incumbent_unique_peptide_conflict')
    if nk and rk and nr>rr+p['risk_advantage']: reasons.append('local_risk_regression')
    unsupported=[iv for iv in ni if not(support.supported_junction(m,iv,votes)
                  or iv in pn['bridged_junctions'])]
    if c.graph_path and unsupported: reasons.append('unsupported_new_junction')
    if delta<=p['replacement_margin']: reasons.append('absolute_evidence_margin')
    return dict(event='replacement',candidate=c.uid,incumbents=[old.uid],
                accepted=not reasons,reasons=reasons,policy=p['selection_policy'],
                channels={k:dict(new=a,old=b,weight=w,delta=w*(a-b)) for k,(a,b,w) in channels.items()},
                normalized_delta=delta,peptide_new=pn,peptide_old=po,
                disputed_intervals=regions,new_junctions_supported=not unsupported,
                triggers=dict(molecular=molecular,residue=residue,rna=rna,protein=protein,
                              peptide=pn['strong'],structural=structural),
                graph_candidate=c.graph_path,missing_evidence_is_negative=False)


def _anchors(originals, source, support, chains, policy):
    """Greedy non-overlapping original anchors, never based on sparse peptides.

    Evidence-dependent changes are tested later by compare_adaptive. Original
    anchors absent from the fallback predictor can be retained in separate
    regions without requiring mass-spec detection.
    """
    index=CodingIndex();events=[]
    # Preserve every fallback original, including genuine overlapping genes.
    # Local alternatives must pass compare_adaptive, rather than bypassing its
    # protections via an evidence-sorted initial anchor.
    for c in sorted(originals,key=lambda c:(source not in c.sources,c.uid)):
        if source in c.sources or (not index.conflicts(c) and c.complete and c.canonical):
            index.add(c)
            events.append(dict(event='original_anchor',candidate=c.uid,accepted=True,
                               reasons=[],policy=policy,missing_evidence_is_negative=False))
    return index,events


def select_adaptive(candidates,source_models,peptides,support,p):
    chains,votes,exons=predictor_support(candidates)
    source,report=choose_backbone(candidates,source_models,chains,p)
    if p['selection_policy']=='adaptive_locus':
        # A genome/consensus fallback holds source choice independent of which
        # experimental evidence is supplied. Evidence acts in local disputes.
        # Fractions, not raw gene/peptide counts, avoid rewarding larger sets.
        originals_by_source={s:[c for c in candidates if s in c.sources and not c.graph_path]
                             for s in source_models}
        source=min(source_models,key=lambda s:(
            -sum(c.complete and c.canonical for c in originals_by_source[s])/max(1,len(originals_by_source[s])),
            -report['sources'][s]['mean_chain_agreement'],s))
    report.update(policy_used=p['selection_policy'],fallback_source=source,
                  anchor_policy='source_then_private_originals' if p['selection_policy']=='adaptive_guarded' else 'genome_consensus_fallback_then_local_evidence',
                  peptide_counts_choose_backbone=False)
    originals=[c for c in candidates if c.sources and not c.graph_path]
    index,events=_anchors(originals,source,support,chains,p['selection_policy'])
    ordered=sorted(candidates,key=lambda c:(c.graph_path,-int(c.complete and c.canonical),
                   -(c.q if c.q is not None else -1),-len(chains.get(chain_key(c.model),())),c.uid))
    for c in ordered:
        if c.uid in index.items: continue
        olds=index.conflicts(c)
        if len(olds)>1:
            events.append(dict(event='replacement',candidate=c.uid,incumbents=[x.uid for x in olds],accepted=False,reasons=['fusion_or_multiple_incumbents']))
        elif olds:
            e=compare_adaptive(c,olds[0],peptides,support,chains,votes,len(source_models),p)
            events.append(e)
            if e['accepted']:index.remove(olds[0]);index.add(c)
        else:
            # Novel models need explicit support. Nonoverlapping original models
            # are retained as anchors, even if no peptide was observed.
            new_ok=all(support.supported_junction(c.model,iv,votes) for iv in c.model.introns)
            supported=(len(chains.get(chain_key(c.model),()))>=2 or
                       (bool(c.model.introns) and support.rna_score(c.model)>=.5) or
                       len(c.peptides)>=2 or support.protein_score(c.model)>=.8)
            accepted=c.complete and c.canonical and new_ok and supported
            events.append(dict(event='addition',candidate=c.uid,incumbents=[],accepted=accepted,
                               reasons=[] if accepted else ['unsupported_novel_addition']))
            if accepted:index.add(c)
    return sorted(index.items.values(),key=lambda c:c.uid),events,source_report(report)


def source_report(report):
    report['experimental']=True
    report['reference_used']=False
    return report
