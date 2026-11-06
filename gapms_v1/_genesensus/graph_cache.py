"""Content-addressed, bounded graph geometry cache independent of score channels."""

from pathlib import Path
from dataclasses import asdict
from concurrent.futures import ProcessPoolExecutor
import gzip, json, os, uuid
from .genome import Genome
from .graph import loci, reconstruct
from .models import Transcript
from .provenance import digest_object, sha256, verify_inputs
from .config import bam_index_paths

_GENOME = None


def _initialize(path, index):
    global _GENOME
    _GENOME = Genome.from_index(path, index)


def _work(args):
    models, junctions, segments, settings, molecular, observations = args
    base,report=reconstruct(models, _GENOME, junctions, extra_segments=segments, molecular_corrections=molecular, **settings)
    if observations:
        guided,guided_report=reconstruct(models,_GENOME,junctions,extra_segments=segments,
            molecular_corrections=molecular,peptide_observations=observations,**settings)
        # Keep every baseline path; peptide guidance cannot erase its search space.
        bykey={m.key:m for m in base}
        bykey.update((m.key,m) for m in guided)
        report['peptide_guided_search']=guided_report
        report['baseline_candidates_preserved']=len(base)
        base=sorted(bykey.values(),key=lambda m:m.key)
        report['candidates']=len(base)
    return base,report


def _decode(row):
    d = dict(row)
    d["cds"] = tuple(tuple(iv) for iv in d["cds"])
    d["exons"] = tuple(tuple(iv) for iv in d["exons"])
    return Transcript(**d)


def graph_candidates(pool, config, fingerprints, support, genome, cache):
    p = config["parameters"]
    settings = dict(
        beam=p["graph_beam"],
        max_candidates=p["graph_candidates"],
        max_nodes=p["graph_max_nodes"],
        max_starts=p["graph_max_starts"],
    )
    originals = {
        source: dict(
            annotation=fingerprints[spec["annotation"]],
            phase_sequence=(
                fingerprints[spec["protein_fasta"]]
                if spec.get("protein_fasta")
                else (
                    fingerprints[spec["residue_profiles"]]
                    if spec.get("residue_profiles")
                    else None
                )
            ),
        )
        for source, spec in config["predictions"].items()
    }
    rna = config.get("rna", {})
    rna_key = {
        k: (fingerprints[v] if k in ("bam", "junctions") else v) for k, v in rna.items()
    }
    if "bam" in rna:
        rna_key["indexes"] = {
            path: fingerprints[path] for path in bam_index_paths(rna["bam"])
        }
    identity = digest_object(
        dict(
            genome=fingerprints[config["genome"]],
            originals=originals,
            rna=rna_key,
            peptides=({k:fingerprints[v] for k,v in config.get('peptides',{}).items()}
                if p.get('peptide_strategy')=='graph_contrast' else None),
            minimum_rna_templates=p["minimum_rna_templates"],
            protein=(
                fingerprints[config["proteins"]["annotation"]]
                if "proteins" in config
                else None
            ),
            settings=settings,
            molecular_guidance=p["selection_policy"]=="state_molecular_guided",
            selection_policy=p["selection_policy"],
            molecular={s:{k:fingerprints[v] for k,v in spec.items() if k in ("scores","residue_profiles")} for s,spec in config["predictions"].items()},
            molecular_parameters={k:p[k] for k in ("localizer_gate","risk_strength","boundary_window_nt","molecular_weight")},
            code={
                name: sha256(Path(__file__).with_name(name + ".py"))
                for name in [
                    "graph_cache",
                    "local_molecular",
                    "molecular",
                    "graph",
                    "genome",
                    "annotation",
                    "normalization",
                    "models",
                    "support",
                    "rna",
                    "peptide_context",
                    "peptides",
                    "provenance",
                ]
            },
        )
    )
    cachepath = Path(cache) / ("graph-" + identity + ".json.gz")
    if cachepath.is_file():
        with gzip.open(cachepath, "rt") as f:
            data = json.load(f)
        if (
            data["identity"] != identity
            or digest_object(data["payload"]) != data["payload_digest"]
        ):
            raise ValueError("Graph cache integrity check failed: " + str(cachepath))
        novel = [(_decode(x["model"]), x["parents"]) for x in data["payload"]["novel"]]
        return novel, data["payload"]["reports"], dict(reused=True, key=identity)
    components = list(loci([c.model for c in pool.values()]))
    from .local_molecular import LocalMolecular
    bykey={c.model.key:c for c in pool.values()}
    args = [
        (mm, support.graph_junctions(mm), support.graph_segments(mm), settings,
         LocalMolecular([bykey[m.key] for m in mm],p).graph_corrections(mm)
         if p['selection_policy']=='state_molecular_guided' and any(bykey[m.key].q is not None for m in mm) else {},
         support.peptide_evidence.graph_observations(mm,support)
         if p.get('peptide_strategy')=='graph_contrast' and hasattr(support,'peptide_evidence') else ())
        for mm in components
    ]
    novel = []
    reports = []

    def collect(results):
        for mm, (found, report) in zip(components, results):
            reports.append(
                dict(
                    contig=mm[0].contig,
                    strand=mm[0].strand,
                    start=min(m.cds[0][0] for m in mm),
                    **report,
                )
            )
            parents = sorted({m.uid for m in mm})
            novel.extend((m, parents) for m in found)

    if p["threads"] == 1:
        _initialize(config['genome'],genome.index)
        collect(_work(row) for row in args)
    else:
        with ProcessPoolExecutor(
            max_workers=p["threads"],
            initializer=_initialize,
            initargs=(config["genome"], genome.index),
        ) as executor:
            collect(executor.map(_work, args, chunksize=16))
    verify_inputs(fingerprints, "graph cache")
    payload = dict(
        novel=[dict(model=asdict(m), parents=parents) for m, parents in novel],
        reports=reports,
    )
    data = dict(
        identity=identity, payload=payload, payload_digest=digest_object(payload)
    )
    temp = cachepath.with_name(
        "." + cachepath.name + ".tmp." + str(os.getpid()) + "." + uuid.uuid4().hex
    )
    with gzip.open(temp, "wt", compresslevel=1) as f:
        json.dump(data, f, sort_keys=True, separators=(",", ":"), allow_nan=False)
    os.replace(temp, cachepath)
    return novel, reports, dict(reused=False, key=identity)
