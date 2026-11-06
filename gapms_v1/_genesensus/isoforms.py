"""Add source-supplied secondary coding isoforms after immutable primary selection."""
from collections import defaultdict, Counter
from .models import Candidate
from .intervals import CodingIndex
from .molecular import read_scores
from .isoform_guards import relationship, routes, qualifies
from .annotation import write_gff3, attributes


def admit_isoforms(primary, source_models, specs, support, genome, parameters):
    """Fixed event_guarded rule; optional input names do not define eligibility."""
    if parameters['isoform_policy'] == 'off':
        return [], {}, [], {'policy': 'off', 'added': 0}
    scores_present = any(spec.get('scores') for spec in specs.values())
    rna_present = parameters.get('_isoform_rna_available', bool(getattr(support, 'rna', None)))
    if not scores_present and not rna_present:
        return [], {}, [], dict(policy='event_guarded', primary_models=len(primary), added=0, primary_geometries_preserved=True, RNA_available=False, scores_available=False, reason='no_isoform_admission_evidence')
    eligible = set()
    for source, models in source_models.items():
        groups = defaultdict(list)
        for m in models:
            groups[m.contig, m.strand, m.gene].append(m)
        for members in groups.values():
            if len({m.cds for m in members}) > 1:
                eligible.update(m.key for m in members)
    pool = {}
    for source, models in source_models.items():
        scores = read_scores(specs[source]['scores']) if specs[source].get('scores') else {}
        for m in models:
            if m.key not in eligible:
                continue
            c = pool.setdefault(m.key, Candidate(m))
            c.sources.add(source)
            c.original_ids.setdefault(source, []).append(m.identifier)
            if m.identifier in scores:
                c.gene_scores.append(scores[m.identifier])
    use_scores = any(spec.get('scores') for spec in specs.values())
    use_rna = bool(getattr(support, 'rna', None))
    # Availability is passed explicitly by the pipeline, even for empty RNA files.
    use_rna = parameters.get('_isoform_rna_available', use_rna)
    index = CodingIndex()
    for c in primary:
        index.add(c)
    selected = {c.uid for c in primary}
    options = []
    rejected = Counter()
    for c in pool.values():
        if c.uid in selected:
            continue
        m = c.model
        complete, canonical, _ = genome.assess(m)
        c.complete, c.canonical = complete, canonical
        if not (complete and canonical and m.length >= 150 and all(b-a >= 20 for a,b in m.introns)):
            rejected['incomplete_noncanonical_or_short'] += 1
            continue
        conflicts = index.conflicts(c)
        if len(conflicts) != 1:
            rejected['zero_or_multiple_primary_genes'] += 1
            continue
        parent = conflicts[0].model
        ok, reason, fraction = relationship(m, parent)
        if not ok:
            rejected[reason] += 1
            continue
        f = routes(c, parent, support, use_scores, use_rna, [], gate=parameters['isoform_score_gate'])
        # This frozen addon does not use a residue-risk veto. Primary selection does.
        rank = (int(f['RNA_event_supported']), min(f['alternative_junction_templates'], default=0) if use_rna else 0, f['q'] if use_scores and f['q'] is not None else -1, len(c.sources), fraction)
        options.append((rank, c, parent, f))
    added, parents, records = [], {}, []
    perparent = defaultdict(list)
    for rank, c, parent, f in sorted(options, key=lambda x:(tuple(-v for v in x[0]), x[1].uid)):
        accepted = qualifies(f, 'event_guarded')
        reason = 'supported' if accepted else 'insufficient_isoform_specific_evidence'
        if accepted and len(perparent[parent.uid]) >= 2:
            accepted, reason = False, 'predeclared_two_additional_isoform_budget'
        if accepted and any(c.model.introns == m.introns for m in perparent[parent.uid]):
            accepted, reason = False, 'duplicate_alternative_intron_chain'
        if accepted:
            perparent[parent.uid].append(c.model)
            added.append(c)
            parents[c.uid] = parent.uid
        records.append(dict(event='isoform', candidate=c.uid, parent=parent.uid, sources=sorted(c.sources), accepted=accepted, reason=reason, **f))
    return added, parents, records, dict(policy='event_guarded', primary_models=len(primary), added=len(added), genes_with_additions=len(perparent), structurally_eligible_candidates=len(options), structural_rejections=dict(rejected), primary_geometries_preserved=True, RNA_available=use_rna, scores_available=use_scores, peptide_admission_route=False)


def write_isoform_gff3(models, path, parents):
    """Keep stable transcript IDs, share primary gene parent, and expand its span."""
    available = {m.uid for m in models}
    parents = {uid: parent for uid, parent in parents.items() if uid in available and parent in available}
    write_gff3(models, path)
    if not parents:
        return
    lookup = {m.uid:m for m in models}
    spans = {}
    for uid, parent in parents.items():
        p, m = lookup[parent], lookup[uid]
        ee = (p.exons or p.cds) + (m.exons or m.cds) + p.cds + m.cds
        lo, hi = spans.get(parent, (min(a for a,b in ee)+1, max(b for a,b in ee)))
        spans[parent] = (min(lo, min(a for a,b in ee)+1), max(hi, max(b for a,b in ee)))
    lines = []
    for line in path.read_text().splitlines():
        if line.startswith('#'):
            lines.append(line)
            continue
        f = line.split('\t'); a = attributes(f[8])
        if f[2] == 'gene':
            uid = a['ID'][4:]
            if uid in parents:
                continue
            if uid in spans:
                f[3], f[4] = map(str, spans[uid])
        elif f[2] == 'mRNA':
            uid = a['Parent'][4:]
            if uid in parents:
                f[8] = f[8].replace(';Parent=GS.g'+uid+';', ';Parent=GS.g'+parents[uid]+';')
        lines.append('\t'.join(f))
    path.write_text('\n'.join(lines)+'\n')
