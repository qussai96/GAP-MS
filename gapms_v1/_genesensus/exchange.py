"""Local, source-symmetric evidence comparisons; no reference features or fitting.

Predictor votes describe agreement, not independent biological observations or
calibrated probabilities. Coding and junction agreement are separate from q.
"""
from bisect import bisect_left, bisect_right
from functools import partial
import numpy as np
from .intervals import merge, subtract, disputed, intersection_length
from .selection import predictor_support, _peptide_pair, _rna_pair
from .molecular import local_risk
from .scoring import weighted_scores


class LocalVotes:
    def __init__(self, source_models):
        self.n = len(source_models)
        grouped = {}
        self.junctions = {}
        for source, models in source_models.items():
            for m in models:
                grouped.setdefault((source, m.contig, m.strand), []).extend(m.cds)
                for iv in m.introns:
                    self.junctions.setdefault((m.contig, m.strand, *iv), set()).add(source)
        self.coverage = {}
        for key, blocks in grouped.items():
            blocks = merge(blocks)
            self.coverage[key] = blocks, [a for a,b in blocks], [b for a,b in blocks]
        self.sources = tuple(sorted(source_models))

    def coding(self, new, old):
        """Expected coding-state agreement on the same symmetric-difference bases."""
        regions = merge(subtract(new.cds, old.cds) + subtract(old.cds, new.cds))
        if not regions:
            return None, None, False
        scores = [0.0, 0.0]; total = 0; certificate = True
        for a,b in regions:
            bounds = {a,b}; source_blocks = []
            for source in self.sources:
                blocks, starts, ends = self.coverage.get((source,new.contig,new.strand),((),(),()))
                picked = blocks[bisect_right(ends,a):bisect_left(starts,b)]
                source_blocks.append(picked)
                for x,y in picked:bounds.update((max(a,x),min(b,y)))
            points = sorted(bounds)
            for x,y in zip(points,points[1:]):
                if x==y:continue
                v = sum(bool(intersection_length(((x,y),), blocks)) for blocks in source_blocks)/self.n
                on = bool(intersection_length(((x,y),), new.cds))
                sn = v if on else 1-v
                scores[0] += (y-x)*sn; scores[1] += (y-x)*(1-sn)
                total += y-x; certificate &= sn > .5
        return scores[0]/total, scores[1]/total, certificate

    def structural(self, new, old):
        a,b,cert = self.coding(new,old); values=[]; certificates=[]
        if a is not None: values.append((a,b));certificates.append(cert)
        ni,oi=set(new.introns),set(old.introns);ivs=ni^oi
        if ivs:
            ns=[]; os=[]; jc=True
            for iv in sorted(ivs):
                v=len(self.junctions.get((new.contig,new.strand,*iv),()))/self.n
                sn=v if iv in ni else 1-v
                ns.append(sn);os.append(1-sn);jc &= sn>.5
            values.append((sum(ns)/len(ns),sum(os)/len(os)));certificates.append(jc)
        if self.n<2 or not values:return 0.,0.,False,False
        return (sum(x for x,y in values)/len(values),sum(y for x,y in values)/len(values),True,all(certificates))


def compare_exchange(c, old, peptides, support, chains, ivotes, n_sources, p, votes):
    m,o=c.model,old.model; reasons=[]
    regions=disputed(m,o,p['boundary_window_nt'])
    sn,so,sk,certificate=votes.structural(m,o)
    nr,nk=local_risk(c,regions,p['localizer_gate']);orr,ok=local_risk(old,regions,p['localizer_gate'])
    qknown=c.q is not None and old.q is not None
    if qknown:
        mn=c.q*(1-p['risk_strength']*(nr if nk else 0));mo=old.q*(1-p['risk_strength']*(orr if ok else 0))
    elif nk and ok:mn,mo=1-p['risk_strength']*nr,1-p['risk_strength']*orr
    else:mn=mo=0.
    pn,po,newpeps,oldpeps=_peptide_pair(c,old,peptides,regions,p['peptide_saturation'])
    rn,ro=_rna_pair(c,old,support)
    cr=tuple((max(a,x),min(b,y)) for a,b in m.cds for x,y in regions if max(a,x)<min(b,y))
    oc=tuple((max(a,x),min(b,y)) for a,b in o.cds for x,y in regions if max(a,x)<min(b,y))
    an,ao=support.protein_score(m,cr),support.protein_score(o,oc)
    new_junctions=set(m.introns)-set(o.introns)
    supported=all(support.supported_junction(m,iv,ivotes) for iv in new_junctions)
    if not(c.complete and c.canonical):reasons.append('incomplete_or_noncanonical')
    if c.graph_path and not supported:reasons.append('unsupported_new_junction')
    if m.length < p['minimum_retained_fraction']*o.length and not(certificate and supported):reasons.append('uncertified_truncation')
    if oldpeps:reasons.append('incumbent_distinguishing_peptide')
    if ro>0 and rn==0:reasons.append('incumbent_distinguishing_rna')
    if nk and ok and nr>orr+p['risk_advantage']:reasons.append('local_risk_regression')
    qadv=qknown and c.q-old.q>=p['q_advantage']
    radv=nk and ok and orr-nr>=p['risk_advantage']
    advantage=(sk and sn>so) or qadv or (radv and not c.graph_path) or pn>po or rn>ro or an>=ao+.1
    if not advantage:reasons.append('insufficient_distinguishing_evidence')
    active=[qknown or(nk and ok),bool(newpeps or oldpeps),abs(rn-ro)>1e-12,abs(an-ao)>1e-12,sk and abs(sn-so)>1e-12]
    weights=[p['molecular_weight'],p['peptide_weight'],p['rna_weight'],p['protein_weight'],p['predictor_weight']]
    values=np.array([[mn,pn,rn,an,sn],[mo,po,ro,ao,so]],dtype=np.float64)
    scores=weighted_scores(values,[active,active],weights,backend='cpu');delta=float(scores[0]-scores[1])
    if delta<=p['replacement_margin']+1e-12:reasons.append('absolute_margin')
    return dict(candidate=c.uid,incumbents=[old.uid],accepted=not reasons,reasons=reasons,
        comparison_policy='local_exchange',graph_candidate=c.graph_path,disputed_intervals=regions,
        channel_order=['molecular','peptide','rna','protein','local_predictor_agreement'],active_channels=active,
        weights=weights,new_values=values[0].tolist(),old_values=values[1].tolist(),
        new_score=float(scores[0]),old_score=float(scores[1]),absolute_margin=delta,
        molecular_measured_both=qknown,residue_known_both=nk and ok,new_disputed_risk=nr,old_disputed_risk=orr,
        new_unique_peptides=newpeps,old_unique_peptides=oldpeps,new_junctions_supported=supported,
        strict_local_majority_certificate=certificate,
        truncation_exception=bool(m.length<p['minimum_retained_fraction']*o.length and certificate and supported))


def select_exchange(candidates, source_models, peptides, support, p):
    from .selection import select
    return select(candidates,source_models,peptides,support,p,
                  _comparison=partial(compare_exchange,votes=LocalVotes(source_models)))
