"""Deterministic RF2-style conservative additions and local guarded replacements."""

from collections import defaultdict
import statistics
import numpy as np
from .intervals import CodingIndex, disputed, overlap, intersection_length
from .molecular import local_risk
from .scoring import weighted_scores

DEFAULTS = dict(
    isoform_policy="event_guarded",
    isoform_score_gate=0.437721461,
    selection_policy="state_molecular_guided",
    molecular_strategy="source_context",
    peptide_strategy="legacy",
    minimum_coding_nt=150,
    minimum_intron_nt=20,
    predictor_weight=1.0,
    confidence_gate=None,
    compress_diagnostics=False,
    risk_strength=0.25,
    q_advantage=0.1,
    risk_advantage=0.2,
    replacement_margin=0.1,
    minimum_retained_fraction=0.9,
    boundary_window_nt=9,
    localizer_gate=None,
    molecular_weight=1.0,
    peptide_weight=1.0,
    rna_weight=1.0,
    protein_weight=0.5,
    structure_weight=0.25,
    peptide_saturation=2.0,
    addition_q=0.437721461,
    minimum_rna_templates=3,
    graph=True,
    graph_beam=8,
    graph_candidates=32,
    graph_max_nodes=2000,
    graph_max_starts=32,
    threads=1,
    backend="cpu",
)


def chain_key(m):
    return (m.contig, m.strand, m.introns if m.introns else ("single", m.cds))


def predictor_support(candidates):
    chains = defaultdict(set)
    introns = defaultdict(set)
    exons = defaultdict(set)
    for c in candidates:
        m = c.model
        chains[chain_key(m)].update(c.sources)
        for a, b in m.introns:
            introns[m.contig, m.strand, a, b].update(c.sources)
        for a, b in m.cds:
            exons[m.contig, m.strand, a, b].update(c.sources)
    return (
        chains,
        {k: len(v) for k, v in introns.items()},
        {k: len(v) for k, v in exons.items()},
    )


def choose_backbone(candidates, source_models, chains, parameters):
    bykey = {c.model.key: c for c in candidates}
    rank = []
    report = {}
    for source, models in sorted(source_models.items()):
        cc = [bykey[m.key] for m in models]
        quality = [c.q for c in cc if c.q is not None]
        fraction = len(quality) / len(cc)
        agreement = statistics.mean(
            len(chains[chain_key(c.model)]) / len(source_models) for c in cc
        )
        report[source] = dict(
            models=len(models),
            q_coverage=fraction,
            median_q=statistics.median(quality) if quality else None,
            mean_chain_agreement=agreement,
        )
    use_molecular = all(v["q_coverage"] >= 0.95 for v in report.values())
    for source, data in report.items():
        value = data["median_q"] if use_molecular else data["mean_chain_agreement"]
        rank.append((-value, source))
    selected = min(rank)[1]
    return selected, dict(
        policy="median_gene_score" if use_molecular else "mean_chain_agreement",
        sources=report,
        selected=selected,
    )


def _structure(c, n_sources, chains):
    return (
        0.5 * len(chains.get(chain_key(c.model), ())) / max(1, n_sources)
        + 0.25 * int(c.complete)
        + 0.25 * int(c.canonical)
    )


def _peptide_pair(c, old, peptides, regions, saturation):
    new_only = {
        s for s in c.peptides - old.peptides if overlap(peptides[s].blocks, regions)
    }
    old_only = {
        s for s in old.peptides - c.peptides if overlap(peptides[s].blocks, regions)
    }
    return (
        len(new_only) / (len(new_only) + saturation),
        len(old_only) / (len(old_only) + saturation),
        sorted(new_only),
        sorted(old_only),
    )


def _rna_pair(c, old, support):
    ni = set(c.model.introns) - set(old.model.introns)
    oi = set(old.model.introns) - set(c.model.introns)
    return support.rna_score(c.model, ni), support.rna_score(old.model, oi)


def compare(c, old, peptides, support, chains, ivotes, n_sources, p):
    m, o = c.model, old.model
    reasons = []
    if not (c.complete and c.canonical):
        reasons.append("incomplete_or_noncanonical")
    if m.length < p["minimum_retained_fraction"] * o.length:
        reasons.append("truncation")
    regions = disputed(m, o, p["boundary_window_nt"])
    nr, nk = local_risk(c, regions, p["localizer_gate"])
    oldr, ok = local_risk(old, regions, p["localizer_gate"])
    qknown = c.q is not None and old.q is not None
    if qknown:
        mn = c.q * (1 - p["risk_strength"] * (nr if nk else 0))
        mo = old.q * (1 - p["risk_strength"] * (oldr if ok else 0))
    elif nk and ok:
        mn = 1 - p["risk_strength"] * nr
        mo = 1 - p["risk_strength"] * oldr
    else:
        mn = mo = 0.0
    pn, po, newpeps, oldpeps = _peptide_pair(
        c, old, peptides, regions, p["peptide_saturation"]
    )
    rn, ro = _rna_pair(c, old, support)
    # Compare protein evidence only over portions of each CDS in the dispute.
    cr = tuple(
        (max(a, x), min(b, y))
        for a, b in m.cds
        for x, y in regions
        if max(a, x) < min(b, y)
    )
    orr = tuple(
        (max(a, x), min(b, y))
        for a, b in o.cds
        for x, y in regions
        if max(a, x) < min(b, y)
    )
    an, ao = support.protein_score(m, cr), support.protein_score(o, orr)
    sn, so = _structure(c, n_sources, chains), _structure(old, n_sources, chains)
    active = [
        qknown or (nk and ok),
        bool(newpeps or oldpeps),
        max(rn, ro) > 0 and abs(rn - ro) > 1e-12,
        max(an, ao) > 0 and abs(an - ao) > 1e-12,
        abs(sn - so) > 1e-12,
    ]
    weights = [
        p["molecular_weight"],
        p["peptide_weight"],
        p["rna_weight"],
        p["protein_weight"],
        p["structure_weight"],
    ]
    qadv = qknown and c.q - old.q >= p["q_advantage"]
    radv = nk and ok and oldr - nr >= p["risk_advantage"]
    new_introns = set(m.introns) - set(o.introns)
    supported = all(support.supported_junction(m, iv, ivotes) for iv in new_introns)
    if c.graph_path and not supported:
        reasons.append("unsupported_new_junction")
    structural_advantage = (
        c.complete
        and not old.complete
        and supported
        and len(chains.get(chain_key(m), ())) >= 2
    )
    advantage = (
        qadv
        or (radv and not c.graph_path)
        or pn > po
        or rn > ro
        or an >= ao + 0.1
        or structural_advantage
    )
    if not advantage:
        reasons.append("insufficient_distinguishing_evidence")
    if nk and ok and nr > oldr + p["risk_advantage"]:
        reasons.append("local_risk_regression")
    values = np.array([[mn, pn, rn, an, sn], [mo, po, ro, ao, so]], dtype=np.float64)
    scores = weighted_scores(values, [active, active], weights, backend="cpu")
    margin = (scores[0] - scores[1]) / max(scores[1], 1e-9)
    if margin <= p["replacement_margin"] + 1e-12:
        reasons.append("margin")
    return dict(
        candidate=c.uid,
        incumbents=[old.uid],
        accepted=not reasons,
        reasons=reasons,
        molecular_measured_both=qknown,
        residue_known_both=nk and ok,
        q_gain=c.q - old.q if qknown else None,
        new_disputed_risk=nr,
        old_disputed_risk=oldr,
        new_unique_peptides=newpeps,
        old_unique_peptides=oldpeps,
        active_channels=active,
        channel_order=["molecular", "peptide", "rna", "protein", "structure"],
        weights=weights,
        new_values=values[0].tolist(),
        old_values=values[1].tolist(),
        weighted_delta=[
            float(w * (a - b)) if on else 0.0
            for a, b, w, on in zip(values[0], values[1], weights, active)
        ],
        new_score=float(scores[0]),
        old_score=float(scores[1]),
        relative_margin=float(margin),
        graph_candidate=c.graph_path,
        disputed_intervals=regions,
        new_junctions_supported=supported,
    )


def select(candidates, source_models, peptides, support, p, _comparison=None, molecular_factory=None):
    if p.get("selection_policy") in ("locus_coverage", "locus_gene", "state_signed", "state_ratio", "state_rna_counts", "state_rna_log", "state_gap", "state_gap_soft", "state_admitted", "state_coherent", "state_soft_admitted", "state_splice_validated", "state_rna_corroborated", "state_molecular_guided") and _comparison is None:
        from .locus_select import select_locus
        return select_locus(candidates, source_models, peptides, support, p, molecular_factory=molecular_factory)
    if p.get("selection_policy") in ("available_exchange", "available_tiered") and _comparison is None:
        from .available import select_available
        return select_available(candidates, source_models, peptides, support, p)
    if p.get("selection_policy") in ("joint_exchange", "joint_tiered") and _comparison is None:
        from .joint import select_joint
        return select_joint(candidates, source_models, peptides, support, p)
    if p.get("selection_policy") in ("molecular_exchange", "molecular_tiered") and _comparison is None:
        from .molecular_exchange import select_molecular
        return select_molecular(candidates, source_models, peptides, support, p)
    if p.get("selection_policy") in ("corroborated_exchange", "corroborated_tiered") and _comparison is None:
        from .corroborated import select_corroborated
        return select_corroborated(candidates, source_models, peptides, support, p)
    if p.get("selection_policy") == "local_exchange" and _comparison is None:
        from .exchange import select_exchange
        return select_exchange(candidates, source_models, peptides, support, p)
    if p.get("selection_policy", "legacy") == "scoped" and _comparison is None:
        from .scoped import compare_scoped
        return select(candidates, source_models, peptides, support, p, _comparison=compare_scoped)
    if p.get("selection_policy", "legacy") != "legacy" and _comparison is None:
        from .adaptive import select_adaptive
        return select_adaptive(candidates, source_models, peptides, support, p)
    chains, ivotes, evotes = predictor_support(candidates)
    source, source_report = choose_backbone(candidates, source_models, chains, p)
    lookup = {c.model.key: c for c in candidates}
    index = CodingIndex()
    events = []
    for m in source_models[source]:
        index.add(lookup[m.key])
    ordered = sorted(
        candidates,
        key=lambda c: (
            c.q is None,
            -(c.q if c.q is not None else 0.0),
            -c.graph_evidence.get("global_evidence_score", 0.0),
            -_structure(c, len(source_models), chains),
            c.uid,
        ),
    )
    for c in ordered:
        if c.uid in index.items:
            continue
        olds = index.conflicts(c)
        m = c.model
        if not olds:
            reasons = []
            if not (c.complete and c.canonical):
                reasons.append("incomplete_or_noncanonical")
            consensus = len(chains.get(chain_key(m), ())) >= 2
            empirical = (
                len(c.peptides) >= 2
                or (bool(m.introns) and support.rna_score(m) >= 0.5)
                or support.protein_score(m) >= 0.8
            )
            independent_exons = any(
                evotes.get((m.contig, m.strand, a, b), 0) >= 2 for a, b in m.cds
            )
            molecular_private = (
                c.q is not None and c.q >= p["addition_q"] and independent_exons
            )
            if not (consensus or empirical or molecular_private):
                reasons.append("unsupported_addition")
            e = dict(
                event="addition",
                candidate=c.uid,
                incumbents=[],
                accepted=not reasons,
                reasons=reasons,
                consensus=consensus,
                empirical=empirical,
                molecular_private=molecular_private,
                graph_candidate=c.graph_path,
            )
            events.append(e)
            if not reasons:
                index.add(c)
        elif len(olds) > 1:
            events.append(
                dict(
                    event="replacement",
                    candidate=c.uid,
                    incumbents=[o.uid for o in olds],
                    accepted=False,
                    reasons=["fusion_or_multiple_incumbents"],
                    graph_candidate=c.graph_path,
                )
            )
        else:
            e = (_comparison or compare)(
                c, olds[0], peptides, support, chains, ivotes, len(source_models), p
            )
            e["event"] = "replacement"
            events.append(e)
            if e["accepted"]:
                index.remove(olds[0])
                index.add(c)
    return sorted(index.items.values(), key=lambda c: c.uid), events, source_report
