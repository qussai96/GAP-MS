"""Public inference API; optional evidence, graph reconstruction and guarded selection."""

from __future__ import annotations
from pathlib import Path
from dataclasses import dataclass
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import time, json, csv, os, resource, traceback, uuid, fcntl, logging
import numpy as np
from . import __version__
from .config import load_config, input_paths, bam_index_paths
from .provenance import (
    sha256,
    digest_object,
    atomic_json,
    software_record,
    verify_inputs,
)
from .annotation import read_annotation, write_gff3
from .genome import Genome
from .models import Candidate
from .diagnostics import diagnostic_writer
from .normalization import normalize_phases
from .loading import load_sources
from .molecular import attach_molecular, local_risk, RiskTrack
from .peptides import read_bed, attach_peptides
from .support import Support
from .graph import loci, reconstruct
from .graph_cache import graph_candidates
from .selection import select, predictor_support, _structure
from .scoring import weighted_scores
from .intervals import intersection_length
from .isoforms import admit_isoforms, write_isoform_gff3


@dataclass(frozen=True)
class RunResult:
    output: Path
    annotation: Path
    report: dict
    reused: bool = False


LOGGER = logging.getLogger("genesensus")
_GRAPH_GENOME = None


def _graph_initialize(path, index):
    global _GRAPH_GENOME
    _GRAPH_GENOME = Genome.from_index(path, index)


def _graph_work(args):
    models, junctions, settings = args
    return reconstruct(models, _GRAPH_GENOME, junctions, **settings)


def _frame_at(m, pos):
    used = 0
    for a, b in m.cds:
        if a <= pos < b:
            offset = (
                used + pos - a if m.strand == "+" else m.length - used - (pos - a) - 1
            )
            return (offset - m.phase) % 3
        used += b - a
    return None


@dataclass
class InheritedRisk:
    track: object
    intervals: tuple

    def risk(self, regions, gate=None):
        intersections = tuple(
            (max(a, x), min(b, y))
            for a, b in self.intervals
            for x, y in regions
            if max(a, x) < min(b, y)
        )
        return self.track.risk(intersections, gate)


def inherit_local_risks(candidate, original_candidates):
    """Same-frame sequence priors only; never manufacture a novel protein q."""
    m = candidate.model
    for original in original_candidates:
        o = original.model
        if (m.contig, m.strand) != (o.contig, o.strand):
            continue
        valid = []
        for a, b in m.cds:
            for x, y in o.cds:
                left, right = max(a, x), min(b, y)
                if left < right and _frame_at(m, left) == _frame_at(o, left):
                    valid.append((left, right))
        if valid:
            for track in original.risks:
                candidate.risks.append(InheritedRisk(track, tuple(valid)))


def _rna_data(config, cache, threads, genome, fingerprints):
    if "rna" not in config:
        return None, dict(available=False)
    from .rna import extract_junctions, write_junctions

    geometry_mode=config.get("parameters",{}).get("selection_policy","legacy") in ("state_soft_admitted","state_splice_validated", "state_rna_corroborated", "state_molecular_guided")
    spec = dict(config["rna"])
    library = spec.pop("library", "unstranded")
    if "junctions" in spec:
        data = json.loads(Path(spec["junctions"]).read_text())
        if data.get("format_version") != 1 or not data.get("references"):
            raise ValueError(
                "RNA junction cache requires format_version=1 and reference lengths"
            )
        if data.get("settings", {}).get("library", library) != library:
            raise ValueError("RNA cache library mode differs from requested library")
        for name, length in data.get("references", {}).items():
            if genome.lengths.get(name) != length:
                raise ValueError("Cached RNA/genome reference mismatch: " + name)
        for row in data.get("junctions", []):
            if (
                len(row) != 5
                or row[0] not in genome.lengths
                or any(isinstance(v, bool) or not isinstance(v, int) for v in row[1:3])
                or not 0 <= row[1] < row[2] <= genome.lengths[row[0]]
            ):
                raise ValueError("Cached RNA junction/genome coordinate mismatch")
        if geometry_mode and data.get('junctions') and not all(k in data for k in ['rna_chain_junctions','observed_internal_exon_counts']):
            raise ValueError('This state policy requires observed RNA-chain geometry: supply an indexed BAM or a geometry cache produced by this release')
        declared_minimum=data.get('alignment_geometry_provenance',{}).get('internal_exon_min_templates')
        if geometry_mode and declared_minimum is not None and declared_minimum!=config['parameters']['minimum_rna_templates']:
            raise ValueError('RNA cache discovery template threshold differs from the run configuration; regenerate the cache with matching settings')
        if config['parameters']['selection_policy'] in ('state_rna_corroborated','state_molecular_guided') and data.get('junctions') and 'observed_internal_exon_counts_all' not in data:
            raise ValueError('RNA corroboration requires unfiltered observed internal-exon counts; use an indexed BAM or a cache produced by this release')
        return data, dict(
            available=True, cached=True, records=len(data.get("junctions", []))
        )
    bam = spec.pop("bam")
    expected = {path: fingerprints[path] for path in [bam] + bam_index_paths(bam)}
    geometry_identity={}
    if geometry_mode:
        expected[str(genome.path)]=fingerprints[str(genome.path)]
        geometry_identity=dict(helper=sha256(Path(__file__).parent/'rna_geometry.py'),minimum_templates=config['parameters']['minimum_rna_templates'])
    key = digest_object(
        dict(
            geometry=geometry_identity,
            inputs=expected,
            references=genome.lengths,
            library=library,
            settings=spec,
            extractor=sha256(Path(__file__).parent / "rna.py"),
        )
    )
    path = cache / f"rna-{key}.json"
    if path.is_file():
        record = json.loads(path.read_text())
        if record.get("identity") != key or digest_object(
            record["payload"]
        ) != record.get("payload_digest"):
            raise ValueError("RNA cache integrity check failed: " + str(path))
        data = record["payload"]
        return data, dict(available=True, cached=True, records=len(data["junctions"]))
    data = extract_junctions(
        bam, genome.lengths, threads=threads, library=library, **spec
    )
    if geometry_mode:
        from .rna_geometry import augment_geometry
        data=augment_geometry(bam,genome,data,threads=threads,minimum_templates=config['parameters']['minimum_rna_templates'])
    verify_inputs(expected, "RNA cache")
    atomic_json(
        path, dict(identity=key, payload=data, payload_digest=digest_object(data))
    )
    return data, dict(available=True, cached=False, records=len(data["junctions"]))


def _run(config_path, output, threads=None, cache=None, resume=False):
    start = time.perf_counter()
    LOGGER.info("Validating configuration and fingerprinting inputs")
    config = load_config(config_path)
    if threads is not None:
        if threads < 1:
            raise ValueError("threads must be positive")
        config["parameters"]["threads"] = threads
    p = config["parameters"]
    software = software_record()
    fingerprints = {path: sha256(path) for path in input_paths(config)}
    identity = digest_object(
        dict(
            config=config,
            inputs=fingerprints,
            version=__version__,
            code=software["package_code"],
        )
    )
    output = Path(output).resolve()
    if output.exists():
        frozen = output / "provenance.json"
        if not resume or not frozen.is_file():
            raise ValueError(
                "Output already exists; use --resume only for a completed identical run"
            )
        record = json.loads(frozen.read_text())
        if record["run_identity"] != identity:
            raise ValueError(
                "Cannot resume: configuration, inputs or package code changed"
            )
        for rel, h in record["output_sha256"].items():
            if sha256(output / rel) != h:
                raise ValueError(
                    "Cannot resume: completed output checksum changed: " + rel
                )
        return RunResult(
            output,
            output / "annotation.gff3",
            json.loads((output / "report.json").read_text()),
            True,
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    work = output.with_name(
        "." + output.name + ".partial." + str(os.getpid()) + "." + uuid.uuid4().hex[:10]
    )
    if work.exists():
        raise ValueError("Partial directory already exists: " + str(work))
    work.mkdir()
    cache = Path(cache).resolve() if cache else output.parent / ".genesensus-cache"
    cache.mkdir(parents=True, exist_ok=True)
    timings = {"fingerprint_inputs": time.perf_counter() - start}
    mark = time.perf_counter()
    try:
        with Genome(config["genome"]) as genome:
            source_models, pool, phase_report, molecular_report, source_cache_report = (
                load_sources(config, genome, fingerprints, cache)
            )
            LOGGER.info(
                "Loaded %d original candidates; molecular validation complete",
                len(pool),
            )
            original_candidates = list(pool.values())
            isoform_source_models = {s: list(mm) for s, mm in source_models.items()}
            protein_models = (
                read_annotation(
                    config["proteins"]["annotation"],
                    "protein",
                    origin="protein_alignment",
                )
                if "proteins" in config
                else []
            )
            # Protein alignments supply indexed segments/junctions, not hundreds
            # of redundant, partial homolog paths or terminal boundary candidates.
            protein_count = len(protein_models)
            for model in protein_models:
                if (
                    model.contig not in genome.lengths
                    or model.cds[-1][1] > genome.lengths[model.contig]
                ):
                    raise ValueError(
                        "Protein alignment/genome coordinate mismatch: "
                        + model.identifier
                    )
            for c in pool.values():
                c.complete, c.canonical, _ = genome.assess(c.model)
            timings["load_validate_inputs"] = time.perf_counter() - mark
            mark = time.perf_counter()
            rna, rna_report = _rna_data(
                config, cache, p["threads"], genome, fingerprints
            )
            support = Support(rna, protein_models, p["minimum_rna_templates"])
            peptides = {}
            peptide_report = dict(available=False)
            if 'peptides' in config and p['peptide_strategy'] != 'legacy':
                peptides, peptide_report = read_bed(config['peptides']['bed'],
                    config['peptides'].get('sequences'), keep_ambiguous=True)
                support.peptide_evidence = peptides
                peptide_report['available'] = True
            del protein_models
            timings["rna_evidence"] = time.perf_counter() - mark
            mark = time.perf_counter()
            LOGGER.info(
                "Optional alignment evidence loaded; starting graph reconstruction"
            )
            graph_reports = []
            novel = []
            graph_cache_report = dict(reused=False, disabled=not p["graph"])
            if p["graph"]:
                novel, graph_reports, graph_cache_report = graph_candidates(
                    pool, config, fingerprints, support, genome, cache
                )
                parent_lookup = {c.uid: c for c in pool.values()}
                for model, parent_ids in novel:
                    if model.key in pool:
                        continue
                    c = Candidate(model, complete=True, canonical=True, graph_path=True)
                    originals = [
                        parent_lookup[uid]
                        for uid in parent_ids
                        if parent_lookup[uid].sources
                    ]
                    inherit_local_risks(c, originals)
                    c.graph_evidence = dict(
                        parent_candidates=parent_ids,
                        molecular_score="unmeasured",
                        residue_evidence="same-frame_local_prior_not_a_rescored_protein",
                    )
                    pool[model.key] = c
            LOGGER.info(
                "Graph stage produced %d novel paths; cache reused=%s",
                len(novel),
                graph_cache_report.get("reused"),
            )
            candidates = sorted(pool.values(), key=lambda c: c.uid)
            if p['selection_policy'] in ('state_splice_validated','state_rna_corroborated','state_molecular_guided'):
                from Bio.Seq import Seq
                for candidate in candidates:
                    if candidate.canonical:
                        continue
                    bad=[]
                    for a,b in candidate.model.introns:
                        if b-a<4:
                            bad.append((a,b));continue
                        left=genome.fetch(candidate.model.contig,a,a+2)
                        right=genome.fetch(candidate.model.contig,b-2,b)
                        if candidate.model.strand=='-':
                            left,right=str(Seq(right).reverse_complement()),str(Seq(left).reverse_complement())
                        if (left,right) not in {('GT','AG'),('GC','AG'),('AT','AC')}:
                            bad.append((a,b))
                    candidate.graph_evidence['noncanonical_introns']=bad

            timings["graph_reconstruction"] = time.perf_counter() - mark
            mark = time.perf_counter()
            if "peptides" in config:
                if p['peptide_strategy'] == 'legacy':
                    peptides, peptide_report = read_bed(
                        config["peptides"]["bed"], config["peptides"].get("sequences")
                    )
                    attach_peptides(candidates, peptides, genome)
                else:
                    peptides.attach(candidates, genome)
                peptide_report["available"] = True
            timings["peptide_mapping"] = time.perf_counter() - mark
            mark = time.perf_counter()
            from .molecular_strategy import prepare as prepare_molecular_strategy
            candidates, source_models, molecular_factory, molecular_strategy_report, diagnostic_events, vetoed_orfs = prepare_molecular_strategy(candidates, source_models, config, genome, fingerprints, cache, support)
            if vetoed_orfs:
                write_gff3(vetoed_orfs, work / "provisional_spurious.gff3")
            atomic_json(work / "molecular_diagnostics.json", diagnostic_events)
            # Numeric bulk scoring is the optional CUDA workload; irregular graph and
            # local decisions remain CPU work. Scores do not impute missing values.
            chains, _, _ = predictor_support(candidates)
            features = []
            available = []
            for c in candidates:
                risk, known = local_risk(c, c.model.cds, p["localizer_gate"])
                molecular = (
                    (c.q * (1 - p["risk_strength"] * (risk if known else 0.0)))
                    if c.q is not None
                    else 0.0
                )
                features.append(
                    [
                        molecular,
                        len(c.peptides) / (len(c.peptides) + p["peptide_saturation"]),
                        support.rna_score(c.model),
                        support.protein_score(c.model),
                        _structure(c, len(source_models), chains),
                    ]
                )
                available.append(
                    [
                        c.q is not None,
                        bool(peptides),
                        rna is not None,
                        protein_count > 0,
                        True,
                    ]
                )
            LOGGER.info(
                "Scoring %d candidates and selecting guarded alternatives",
                len(candidates),
            )
            prior = weighted_scores(
                np.asarray(features).reshape(-1,5),
                np.asarray(available).reshape(-1,5),
                [
                    p["molecular_weight"],
                    p["peptide_weight"],
                    p["rna_weight"],
                    p["protein_weight"],
                    p["structure_weight"],
                ],
                p["backend"],
            )
            for c, value in zip(candidates, prior):
                c.graph_evidence["global_evidence_score"] = float(value)
            timings["bulk_scoring"] = time.perf_counter() - mark
            mark = time.perf_counter()
            selected, decisions, backbone = select(
                candidates, source_models, peptides, support, p, molecular_factory=molecular_factory
            )
            timings["guarded_selection"] = time.perf_counter() - mark
            isoform_start = time.perf_counter()
            primary_count = len(selected)
            if p['isoform_policy'] != 'off':
                write_gff3([c.model for c in selected], work / 'primary_annotation.gff3')
            isoform_parameters = dict(p, _isoform_rna_available='rna' in config)
            isoform_added, isoform_parents, isoform_decisions, isoform_report = admit_isoforms(
                selected, isoform_source_models, config['predictions'], support, genome, isoform_parameters)
            selected = list(selected) + isoform_added
            decisions.extend(isoform_decisions)
            known_uids = {c.uid for c in candidates}
            candidates.extend(c for c in isoform_added if c.uid not in known_uids)
            if p['isoform_policy'] != 'off':
                write_gff3([c.model for c in isoform_added], work / 'additional_isoforms.gff3')
                atomic_json(work / 'isoforms.json', dict(summary=isoform_report, decisions=isoform_decisions, parent_by_candidate=isoform_parents))
            timings['isoform_admission'] = time.perf_counter() - isoform_start
            mark = time.perf_counter()
            LOGGER.info(
                "Selected %d models; writing and verifying outputs", len(selected)
            )
            write_isoform_gff3([c.model for c in selected], work / "annotation.gff3", isoform_parents)
            with (work / "proteins.faa").open("w") as f:
                for c in sorted(selected, key=lambda c: c.uid):
                    _, _, protein = genome.assess(c.model)
                    f.write(">GS.g" + c.uid + ".t1\n" + protein.rstrip("*") + "\n")
            deferred = {d['candidate'] for d in decisions if d['event']=='provisional' and d['accepted']}
            deferred -= {c.uid for c in selected}
            if deferred:
                provisional = [c for c in candidates if c.uid in deferred]
                write_gff3([c.model for c in provisional], work / 'provisional.gff3')
                with (work/'provisional_proteins.faa').open('w') as f:
                    for c in sorted(provisional,key=lambda c:c.uid):
                        _,_,protein=genome.assess(c.model)
                        f.write('>GS.g'+c.uid+'.t1\n'+protein.rstrip('*')+'\n')
            if p['selection_policy'] in ('corroborated_tiered','molecular_tiered','joint_tiered','available_tiered'):
                write_isoform_gff3([c.model for c in selected]+[c.model for c in candidates if c.uid in deferred], work/'full_annotation.gff3', isoform_parents)
            if p['selection_policy']=='available_tiered':
                strict = {d['candidate'] for d in decisions if d['event']=='consensus_provisional' and d['accepted']}
                all_selected = list(selected)+[c for c in candidates if c.uid in deferred]
                write_isoform_gff3([c.model for c in all_selected if c.uid not in strict],work/'consensus_annotation.gff3', isoform_parents)
                write_gff3([c.model for c in all_selected if c.uid in strict],work/'consensus_provisional.gff3')
                with (work/'consensus_provisional_proteins.faa').open('w') as f:
                    for c in sorted(all_selected,key=lambda c:c.uid):
                        if c.uid in strict:
                            _,_,protein=genome.assess(c.model)
                            f.write('>GS.g'+c.uid+'.t1\n'+protein.rstrip('*')+'\n')
            ids = {c.uid for c in selected}
            with diagnostic_writer(work / ("candidates.tsv.gz" if p["compress_diagnostics"] else "candidates.tsv"), p["compress_diagnostics"]) as f:
                w = csv.writer(f, delimiter="\t")
                w.writerow(
                    [
                        "candidate",
                        "contig",
                        "strand",
                        "CDS_0based_halfopen",
                        "origin",
                        "sources",
                        "source_ids",
                        "gene_q",
                        "q_measured",
                        "residue_tracks",
                        "complete",
                        "canonical",
                        "peptides",
                        "selected",
                        "graph_details",
                    ]
                )
                for c in candidates:
                    w.writerow(
                        [
                            c.uid,
                            c.model.contig,
                            c.model.strand,
                            json.dumps(c.model.cds),
                            c.model.origin,
                            ",".join(sorted(c.sources)),
                            json.dumps(c.original_ids,sort_keys=True),
                            c.q,
                            c.q is not None,
                            len(c.risks),
                            c.complete,
                            c.canonical,
                            len(c.peptides),
                            c.uid in ids,
                            json.dumps(c.graph_evidence, sort_keys=True),
                        ]
                    )
            with diagnostic_writer(work / ("decisions.jsonl.gz" if p["compress_diagnostics"] else "decisions.jsonl"), p["compress_diagnostics"]) as f:
                for d in decisions:
                    f.write(json.dumps(d, sort_keys=True, allow_nan=False) + "\n")
            atomic_json(work / "graph_report.json", graph_reports)
            timings["write_outputs"] = time.perf_counter() - mark
        mark = time.perf_counter()
        for path, h in fingerprints.items():
            if sha256(path) != h:
                raise ValueError("Input changed while inference was running: " + path)
        for name, h in software["package_code"].items():
            if sha256(Path(__file__).with_name(name)) != h:
                raise ValueError(
                    "Package code changed while inference was running: " + name
                )
        timings["verify_inputs"] = time.perf_counter() - mark
        timings["total_seconds"] = time.perf_counter() - start
        report = dict(
            version=__version__,
            predictors={s: len(mm) for s, mm in source_models.items()},
            candidate_count=len(candidates),
            graph_cache=graph_cache_report,
            selected_count=len(selected),
            primary_selected_count=primary_count,
            isoforms=isoform_report,
            diagnostics_compressed=p["compress_diagnostics"],
            provisional_count=sum(d["event"]=="provisional" and d["accepted"] for d in decisions),
            consensus_provisional_count=sum(d["event"]=="consensus_provisional" and d["accepted"] for d in decisions),
            novel_graph_candidates=sum(c.graph_path for c in candidates),
            selected_graph_candidates=sum(c.graph_path for c in selected),
            backbone=backbone,
            molecular=molecular_report,
            molecular_strategy=molecular_strategy_report,
            source_cache=source_cache_report,
            phase_validation=phase_report,
            peptides=peptide_report,
            rna=rna_report,
            protein_alignment_models=protein_count,
            timings=timings,
            peak_parent_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            max_finished_child_rss_kib=resource.getrusage(
                resource.RUSAGE_CHILDREN
            ).ru_maxrss,
            accepted_additions=sum(
                d["accepted"] and d["event"] == "addition" for d in decisions
            ),
            accepted_replacements=sum(
                d["accepted"] and d["event"] == "replacement" for d in decisions
            ),
            reference_used=False,
            trained=molecular_strategy_report["trained"],
            reference_trained=False,
            backend=p["backend"],
        )
        atomic_json(work / "report.json", report)
        atomic_json(work / "configuration.json", config)
        outputs = {p.name: sha256(p) for p in work.iterdir() if p.is_file()}
        atomic_json(
            work / "provenance.json",
            dict(
                run_identity=identity,
                created_utc=datetime.now(timezone.utc).isoformat(),
                input_sha256=fingerprints,
                output_sha256=outputs,
                software=software,
                reference_used=False,
                trained=molecular_strategy_report["trained"],
            reference_trained=False,
            ),
        )
        if output.exists():
            raise ValueError("Output appeared during run; refusing overwrite")
        os.rename(work, output)
        return RunResult(output, output / "annotation.gff3", report)
    except BaseException as exc:
        atomic_json(
            work / "failure.json",
            dict(
                error=type(exc).__name__,
                message=str(exc),
                traceback=traceback.format_exc(),
                completed=False,
            ),
        )
        raise


def run(config_path, output, threads=None, cache=None, resume=False):
    """Run with an OS-released advisory lock and atomic completed output publication."""
    target = Path(output).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    lockpath = target.parent / ("." + target.name + ".lock")
    with lockpath.open("a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(
                "Another GeneSensus run owns this output: " + str(target)
            ) from exc
        lock.seek(0)
        lock.truncate()
        lock.write(json.dumps(dict(pid=os.getpid(), output=str(target))) + "\n")
        lock.flush()
        try:
            return _run(config_path, target, threads, cache, resume)
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
