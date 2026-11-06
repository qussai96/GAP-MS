"""Parallel, content-addressed source validation and molecular feature caching."""

from pathlib import Path
from dataclasses import asdict
from concurrent.futures import ProcessPoolExecutor
import gzip, json, os, uuid
import numpy as np
from .annotation import read_annotation
from .normalization import normalize_phases
from .molecular import attach_molecular, RiskTrack
from .models import Transcript, Candidate
from .genome import Genome
from .provenance import digest_object, sha256, verify_inputs

_GENOME = None


def _initialize(path, index):
    global _GENOME
    _GENOME = Genome.from_index(path, index)


def _decode(d):
    d = dict(d)
    d["cds"] = tuple(tuple(x) for x in d["cds"])
    d["exons"] = tuple(tuple(x) for x in d["exons"])
    return Transcript(**d)


def _pack(models, candidates, phase, molecular):
    positions = {id(m): i for i, m in enumerate(models)}
    rows = []
    for c in candidates:
        risks = []
        for track in c.risks:
            v = np.asarray(track.probabilities)
            risks.append(
                dict(
                    model_index=positions[id(track.model)],
                    dtype=v.dtype.str,
                    data=v.tobytes().hex(),
                    q=track.q,
                )
            )
        rows.append(
            dict(model_index=positions[id(c.model)], scores=c.gene_scores, risks=risks)
        )
    return dict(
        models=[asdict(m) for m in models],
        candidates=rows,
        phase=phase,
        molecular=molecular,
    )


def _unpack(data):
    models = [_decode(m) for m in data["models"]]
    byuid = {}
    for row in data["candidates"]:
        model = models[row["model_index"]]
        c = Candidate(model)
        byuid[model.uid] = c
        c.gene_scores = list(row["scores"])
        for r in row["risks"]:
            dtype = np.dtype(r["dtype"])
            if dtype.kind != "f" or dtype.itemsize not in (2, 4, 8):
                raise ValueError("Invalid cached probability dtype")
            values = np.frombuffer(bytes.fromhex(r["data"]), dtype=dtype)
            c.risks.append(RiskTrack(models[r["model_index"]], values, r["q"]))
    return models, list(byuid.values()), data["phase"], data["molecular"]


def _source(args, genome):
    source, spec, key, cache, expected = args
    path = Path(cache) / ("source-" + key + ".json.gz")
    if path.is_file():
        with gzip.open(path, "rt") as f:
            data = json.load(f)
        if data["identity"] != key or digest_object(data["payload"]) != data["digest"]:
            raise ValueError("Source cache integrity check failed: " + str(path))
        return source, *_unpack(data["payload"]), dict(reused=True, key=key)
    models = read_annotation(spec["annotation"], source)
    phase = normalize_phases({source: models}, {source: spec}, genome)[source]
    pool = {}
    for m in models:
        pool.setdefault(m.key, Candidate(m))
    molecular = attach_molecular(
        list(pool.values()), {source: models}, {source: spec}, genome
    ).get(source, {})
    verify_inputs(expected, "source cache")
    payload = _pack(models, list(pool.values()), phase, molecular)
    data = dict(identity=key, payload=payload, digest=digest_object(payload))
    temp = path.with_name(
        "." + path.name + ".tmp." + str(os.getpid()) + "." + uuid.uuid4().hex
    )
    with gzip.open(temp, "wt", compresslevel=1) as f:
        json.dump(data, f, sort_keys=True, separators=(",", ":"), allow_nan=False)
    os.replace(temp, path)
    return (
        source,
        models,
        list(pool.values()),
        phase,
        molecular,
        dict(reused=False, key=key),
    )


def _work(args):
    return _source(args, _GENOME)


def load_sources(config, genome, fingerprints, cache):
    code = {
        name: sha256(Path(__file__).with_name(name + ".py"))
        for name in [
            "loading",
            "annotation",
            "normalization",
            "molecular",
            "models",
            "genome",
            "provenance",
        ]
    }
    args = []
    for source, spec in sorted(config["predictions"].items()):
        key = digest_object(
            dict(
                source=source,
                inputs={name: fingerprints[path] for name, path in spec.items()},
                genome=fingerprints[config["genome"]],
                code=code,
            )
        )
        expected = {path: fingerprints[path] for path in spec.values()}
        expected[config["genome"]] = fingerprints[config["genome"]]
        args.append((source, spec, key, str(cache), expected))
    if config["parameters"]["threads"] == 1:
        results = [_source(a, genome) for a in args]
    else:
        with ProcessPoolExecutor(
            max_workers=min(config["parameters"]["threads"], len(args)),
            initializer=_initialize,
            initargs=(config["genome"], genome.index),
        ) as executor:
            results = list(executor.map(_work, args))
    source_models = {}
    pool = {}
    phase_report = {}
    molecular_report = {}
    cache_report = {}
    for source, models, candidates, phase, molecular, cache_info in results:
        source_models[source] = models
        phase_report[source] = phase
        molecular_report[source] = molecular
        cache_report[source] = cache_info
        for m in models:
            c = pool.setdefault(m.key, Candidate(m))
            c.sources.add(source)
            c.original_ids.setdefault(source, []).append(m.identifier)
        for item in candidates:
            c = pool[item.model.key]
            c.gene_scores.extend(item.gene_scores)
            c.risks.extend(item.risks)
    return source_models, pool, phase_report, molecular_report, cache_report
