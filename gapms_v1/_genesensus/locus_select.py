"""Backbone-free local coding-path compatibility optimization; no reference inputs."""
from collections import defaultdict
from functools import lru_cache
import math
from .intervals import overlap,intersection_length
from .graph import loci
from .molecular import local_risk

def similarity(a,b):
    coding=2*intersection_length(a.cds,b.cds)/(a.length+b.length)
    ai,bi=set(a.introns),set(b.introns)
    if not ai and not bi:return coding
    junction=2*len(ai&bi)/(len(ai)+len(bi))
    return (coding+junction)/2

def independent_set(weights,conflicts,exact_limit=28):
    n=len(weights)
    if n<=exact_limit:
        @lru_cache(None)
        def solve(mask):
            if not mask:return (0.,())
            members=[i for i in range(n) if mask>>i&1]
            i=max(members,key=lambda j:((conflicts[j]&mask).bit_count() if hasattr(int,'bit_count') else bin(conflicts[j]&mask).count('1'),weights[j],-j))
            remaining=mask&~(1<<i)
            no=solve(remaining);score,chosen=solve(remaining&~conflicts[i]);yes=(score+weights[i],tuple(sorted(chosen+(i,))))
            return max(no,yes,key=lambda x:(round(x[0],9),tuple(-v for v in x[1])))
        value,chosen=solve((1<<n)-1)
        return list(chosen),dict(exact=True,states=solve.cache_info().currsize,n=n)
    # Bounded beam over independent subsets: honest approximation, retains several orders.
    orders=[sorted(range(n),key=lambda i:(-weights[i],i)),sorted(range(n),key=lambda i:(-weights[i]/(1+bin(conflicts[i]).count('1')),i))]
    answers=[]
    for order in orders:
        chosen=[];mask=0
        for i in order:
            if not conflicts[i]&mask:chosen.append(i);mask|=1<<i
        # Strictly improving single-candidate exchanges; deterministic finite descent.
        changed=True
        while changed:
            changed=False
            for i in order:
                if i in chosen:continue
                old=[j for j in chosen if conflicts[i]>>j&1]
                if weights[i]>sum(weights[j] for j in old)+1e-9:
                    chosen=[j for j in chosen if j not in old]+[i];changed=True
        answers.append((sum(weights[i] for i in chosen),tuple(sorted(chosen))))
    score,chosen=max(answers,key=lambda x:(round(x[0],9),tuple(-v for v in x[1])))
    return list(chosen),dict(exact=False,solver='two deterministic greedy orders with strictly improving exchanges',n=n)

def select_locus(candidates,source_models,peptides,support,p,molecular_factory=None):
    if p["selection_policy"]=="state_molecular_guided":
        from .molecular_state_utility import StateEvidence
    else:
        from .state_utility import StateEvidence
    states=StateEvidence(source_models) if p['selection_policy'].startswith('state_') else None
    from .exon_support import ExonSupport
    exon_support=ExonSupport(source_models,support,peptides) if p['selection_policy'] in ('state_coherent','state_splice_validated','state_rna_corroborated','state_molecular_guided') else None
    bykey={c.model.key:c for c in candidates};selected=[];events=[];summary=defaultdict(int)
    for models in loci([c.model for c in candidates]):
        cc=sorted([bykey[m.key] for m in models],key=lambda c:c.uid);originals=[c for c in cc if c.sources]
        from .local_molecular import LocalMolecular
        molecular_states=LocalMolecular(originals,p) if p["selection_policy"]=="state_molecular_guided" else None
        if molecular_states is not None and p.get("direct_candidate_audit",False):
            from .direct_molecular import DirectMolecular
            molecular_states=DirectMolecular(cc,p)
        if molecular_factory is not None:
            molecular_states=molecular_factory(originals,p)
        weights=[];features=[];valid=[]
        for c in cc:
            m=c.model;votes=[]
            for source in sorted(source_models):
                matches=[similarity(m,o.model) for o in originals if source in o.sources and overlap(m.cds,o.model.cds)]
                votes.append(max(matches,default=0.))
            agreement=sum(votes);risk,known=local_risk(c,m.cds,p['localizer_gate']);molecular=c.q*(1-p['risk_strength']*(risk if known else 0)) if c.q is not None else 0.
            peptide=len(c.peptides)/(len(c.peptides)+p['peptide_saturation']);rna=support.rna_score(m);protein=support.protein_score(m)
            # Coverage objective; unweighted gene-event form is a declared ablation.
            units=m.length if p['selection_policy']=='locus_coverage' else 1000.
            extra=p['molecular_weight']*molecular+p['peptide_weight']*peptide+p['rna_weight']*rna+p['protein_weight']*protein
            score=units*(agreement+extra)
            state_features=(states.features(c,peptides,support,p,molecular_states) if molecular_states is not None else states.features(c,peptides,support,p)) if states else {}
            if states:
                score=state_features['utility']
                if molecular_states is not None:molecular=state_features['molecular']
            # Tie-break quality only, never a material substitute for measured agreement.
            if score>0:score+=1e-7*(int(c.complete)+int(c.canonical))
            reason='noncoding_preferred' if states and state_features['utility']<=0 else None
            if p['selection_policy'] in ('state_admitted','state_coherent','state_soft_admitted','state_splice_validated','state_rna_corroborated','state_molecular_guided') and not state_features['admitted']:
                reason='below_minimum_coding_length' if m.length<p['minimum_coding_nt'] else 'insufficient_coding_to_noncoding_support'
            if p['selection_policy'] in ('state_splice_validated','state_rna_corroborated','state_molecular_guided'):
                if any(b-a<p['minimum_intron_nt'] for a,b in m.introns):
                    reason='intron_below_minimum_length'
                else:
                    if not c.canonical and m.introns and 'noncanonical_introns' not in c.graph_evidence:
                        reason='unvalidated_splice_sites'
                    for iv in c.graph_evidence.get('noncanonical_introns',()):
                        # RNA or a uniquely placed, frame-verified split peptide
                        # can validate an otherwise nonstandard splice site.
                        rna_ok=support.rna_count(m,iv)>=p['minimum_rna_templates']
                        peptide_ok=any(iv in tuple((a[1],b[0]) for a,b in zip(peptides[key].blocks,peptides[key].blocks[1:])) for key in c.peptides)
                        if not (rna_ok or peptide_ok):reason='noncanonical_splice_without_direct_support'
            if c.graph_path:
                if exon_support is not None and not exon_support.supports(m,c.peptides):reason='unsupported_exon_pairing_or_source_frame'
                # Every new junction needs an original predictor or positive empirical support.
                original_introns={iv for o in originals for iv in o.model.introns}
                if any(iv not in original_introns and not support.supported_junction(m,iv) for iv in m.introns):reason='unsupported_novel_junction'
                # A graph path cannot fuse two nonoverlapping multi-source anchors.
                anchors=[o for o in originals if len(o.sources)>=2 and overlap(m.cds,o.model.cds)]
                for i,a in enumerate(anchors):
                    for b in anchors[i+1:]:
                        if a.model.cds[-1][1]<=b.model.cds[0][0] or b.model.cds[-1][1]<=a.model.cds[0][0]:
                            if not all(sum(iv in o.model.introns for o in originals)>=2 or support.rna_count(m,iv)>=p['minimum_rna_templates'] for iv in m.introns):reason='unsupported_multi_anchor_fusion'
            valid.append(reason is None);weights.append(score if reason is None else 0.);features.append(dict(agreement_by_source=votes,molecular=molecular,peptide=peptide,rna=rna,protein=protein,units=units,utility=score,excluded_reason=reason,state_features=state_features))
        edges=[0]*len(cc)
        for i,a in enumerate(cc):
            for j in range(i):
                if overlap(a.model.cds,cc[j].model.cds):edges[i]|=1<<j;edges[j]|=1<<i
        # Solve each actual coding-conflict component; intronic nested genes remain compatible.
        unseen=set(range(len(cc)));chosen=set()
        while unseen:
            todo=[min(unseen)];component=set()
            while todo:
                i=todo.pop()
                if i in component:continue
                component.add(i);unseen.discard(i);todo.extend(j for j in unseen if edges[i]>>j&1)
            indices=sorted(i for i in component if valid[i]);local={g:i for i,g in enumerate(indices)}
            masks=[sum(1<<local[j] for j in indices if edges[i]>>j&1) for i in indices]
            answer,report=independent_set([weights[i] for i in indices],masks)
            chosen.update(indices[i] for i in answer);summary['exact_components' if report['exact'] else 'approximate_components']+=1
        if p.get('peptide_strategy')=='graph_contrast' and hasattr(peptides,'conservative_repair'):
            chosen,peptide_events=peptides.conservative_repair(cc,chosen,valid,weights,support,p)
            events.extend(peptide_events)
            summary['peptide_contrast_repairs']+=len(peptide_events)
        for i,c in enumerate(cc):
            events.append(dict(event='locus_optimization',candidate=c.uid,incumbents=[],accepted=i in chosen,reasons=[] if i in chosen else [features[i]['excluded_reason'] or 'lower_joint_compatibility_utility'],**features[i]))
            if i in chosen:selected.append(c)
            elif c.sources and features[i]['excluded_reason'] is not None:
                events.append(dict(event='provisional',candidate=c.uid,incumbents=[],accepted=True,reasons=[features[i]['excluded_reason']],interpretation='Unresolved source hypothesis retained separately; exclusion is not evidence of gene absence'))
    return selected,events,dict(policy=p['selection_policy'],selected=None,global_backbone=False,sources={},solver_counts=dict(summary))
