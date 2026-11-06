"""Strict inference-only configuration; reference evaluation has a separate interface."""

from pathlib import Path
import json, math, re
from .selection import DEFAULTS

PATH_KEYS = {
    "annotation",
    "scores",
    "residue_profiles",
    "diagnostics",
    "protein_fasta",
    "bed",
    "sequences",
    "bam",
    "junctions",
}


def load_config(path):
    path = Path(path).resolve()
    data = json.loads(path.read_text())
    allowed = {"genome", "predictions", "peptides", "rna", "proteins", "parameters"}
    if not isinstance(data, dict) or set(data) - allowed:
        raise ValueError(
            "Unknown inference configuration keys: " + repr(sorted(set(data) - allowed))
        )
    if (
        not isinstance(data.get("genome"), str)
        or not isinstance(data.get("predictions"), dict)
        or not data["predictions"]
    ):
        raise ValueError("Config requires genome and nonempty named predictions")

    def resolve(p):
        if not isinstance(p, str) or not p:
            raise ValueError("Input paths must be nonempty strings")
        q = Path(p).expanduser()
        q = q if q.is_absolute() else path.parent / q
        if not q.is_file():
            raise ValueError("Input file does not exist: " + str(q))
        return str(q.resolve())

    result = dict(genome=resolve(data["genome"]), predictions={})
    for source, spec in sorted(data["predictions"].items()):
        if not re.fullmatch("[A-Za-z0-9_.-]+", source) or source.lower() in (
            "refseq",
            "reference",
            "araport",
            "truth",
            "ground_truth",
        ):
            raise ValueError("Invalid or reserved predictor source name: " + source)
        if isinstance(spec, str):
            spec = {"annotation": spec}
        if (
            not isinstance(spec, dict)
            or set(spec) - {"annotation", "scores", "residue_profiles", "protein_fasta", "diagnostics"}
            or "annotation" not in spec
        ):
            raise ValueError("Invalid predictor specification: " + source)
        result["predictions"][source] = {k: resolve(v) for k, v in spec.items()}
        for kind, p in result["predictions"][source].items():
            name = Path(p).name.lower()
            if ".ref." in name or "araport" in name or name.startswith("refseq_"):
                raise ValueError(
                    "Target reference annotation/scores are not inference inputs: " + p
                )
    for name, keys, required in [
        ("peptides", {"bed", "sequences"}, "bed"),
        (
            "rna",
            {
                "bam",
                "junctions",
                "library",
                "min_mapq",
                "min_overhang",
                "min_intron",
                "max_intron",
                "skip_duplicates",
            },
            None,
        ),
        ("proteins", {"annotation"}, "annotation"),
    ]:
        spec = data.get(name)
        if spec is None:
            continue
        if (
            not isinstance(spec, dict)
            or set(spec) - keys
            or (required and required not in spec)
        ):
            raise ValueError("Invalid " + name + " evidence specification")
        if name == "rna" and sum(k in spec for k in ("bam", "junctions")) != 1:
            raise ValueError("RNA requires exactly one of bam or junctions")
        result[name] = {k: resolve(v) if k in PATH_KEYS else v for k, v in spec.items()}
    params = dict(DEFAULTS)
    extra = data.get("parameters", {})
    if not isinstance(extra, dict) or set(extra) - set(params):
        raise ValueError(
            "Unknown algorithm parameters: " + repr(sorted(set(extra) - set(params)))
        )
    params.update(extra)
    if params['peptide_strategy'] not in ('legacy','bounded_unique','fractional','graph_contrast'):
        raise ValueError('Invalid peptide_strategy')
    if params['peptide_strategy'] != 'legacy' and params['selection_policy'] != 'state_molecular_guided':
        raise ValueError('Peptide context strategies require state_molecular_guided selection')
    if params["isoform_policy"] not in ("off", "event_guarded"):
        raise ValueError("Invalid isoform_policy")
    if params["molecular_strategy"] not in ("source_context", "splice_context", "splice_orf"):
        raise ValueError("Invalid molecular_strategy")
    if params["molecular_strategy"] != "source_context" and params["selection_policy"] != "state_molecular_guided":
        raise ValueError("Molecular splice strategies require state_molecular_guided selection")
    if params["molecular_strategy"] != "source_context" and params["localizer_gate"] is None:
        params["localizer_gate"] = params["addition_q"]
    if isinstance(params["minimum_intron_nt"],bool) or not isinstance(params["minimum_intron_nt"],int) or params["minimum_intron_nt"]<4:
        raise ValueError("minimum_intron_nt must be an integer >=4")
    if isinstance(params["minimum_coding_nt"],bool) or not isinstance(params["minimum_coding_nt"],int) or params["minimum_coding_nt"]<0:
        raise ValueError("minimum_coding_nt must be a nonnegative integer")
    if params["selection_policy"] not in ("state_signed", "state_ratio", "state_rna_counts", "state_rna_log", "state_gap", "state_gap_soft", "state_admitted", "state_coherent", "state_soft_admitted", "state_splice_validated", "state_rna_corroborated", "state_molecular_guided", "locus_coverage", "locus_gene", "legacy", "adaptive_guarded", "adaptive_locus", "scoped", "local_exchange", "corroborated_exchange", "corroborated_tiered", "molecular_exchange", "molecular_tiered", "joint_exchange", "joint_tiered", "available_exchange", "available_tiered"):
        raise ValueError("Invalid selection_policy")
    for key in [
        "threads",
        "graph_beam",
        "graph_candidates",
        "graph_max_nodes",
        "graph_max_starts",
        "minimum_rna_templates",
    ]:
        if (
            isinstance(params[key], bool)
            or not isinstance(params[key], int)
            or params[key] < 1
        ):
            raise ValueError(key + " must be a positive integer")
    if not isinstance(params["compress_diagnostics"], bool):
        raise ValueError("compress_diagnostics must be boolean")
    if not isinstance(params["graph"], bool):
        raise ValueError("graph must be boolean")
    for key in [
        "isoform_score_gate",
        "risk_strength",
        "q_advantage",
        "risk_advantage",
        "minimum_retained_fraction",
        "addition_q",
    ]:
        if (
            not isinstance(params[key], (int, float))
            or isinstance(params[key], bool)
            or not math.isfinite(params[key])
            or not 0 <= params[key] <= 1
        ):
            raise ValueError(key + " must be in [0,1]")
    if params["minimum_retained_fraction"] == 0:
        raise ValueError("Retained fraction must be positive")
    for key in [
        "replacement_margin",
        "molecular_weight",
        "predictor_weight",
        "peptide_weight",
        "rna_weight",
        "protein_weight",
        "structure_weight",
        "peptide_saturation",
        "boundary_window_nt",
    ]:
        if (
            not isinstance(params[key], (int, float))
            or isinstance(params[key], bool)
            or not math.isfinite(params[key])
            or params[key] < 0
        ):
            raise ValueError(key + " must be finite and nonnegative")
    if params["peptide_saturation"] == 0:
        raise ValueError("Peptide saturation must be positive")
    if params["backend"] not in ("cpu", "auto", "cuda"):
        raise ValueError("Invalid scoring backend")
    if params["confidence_gate"] is not None and (isinstance(params["confidence_gate"], bool) or not isinstance(params["confidence_gate"], (float,int)) or not math.isfinite(params["confidence_gate"]) or not 0 <= params["confidence_gate"] <= 1):
        raise ValueError("confidence_gate must be finite in [0,1]")
    if params["localizer_gate"] is not None and not 0 <= params["localizer_gate"] <= 1:
        raise ValueError("Invalid localizer gate")
    if not isinstance(params["boundary_window_nt"], int) or isinstance(
        params["boundary_window_nt"], bool
    ):
        raise ValueError("boundary_window_nt must be an integer")
    if params["localizer_gate"] is not None and (
        isinstance(params["localizer_gate"], bool)
        or not isinstance(params["localizer_gate"], (int, float))
        or not math.isfinite(params["localizer_gate"])
    ):
        raise ValueError("localizer_gate must be a finite number or null")
    if (
        sum(
            params[k]
            for k in [
                "molecular_weight",
        "predictor_weight",
                "peptide_weight",
                "rna_weight",
                "protein_weight",
                "structure_weight",
            ]
        )
        <= 0
    ):
        raise ValueError("At least one evidence weight must be positive")
    if "rna" in result:
        rna = result["rna"]
        if rna.get("library", "unstranded") not in (
            "unstranded",
            "firststrand",
            "secondstrand",
        ):
            raise ValueError("Invalid RNA library mode")
        for key in ["min_mapq", "min_overhang", "min_intron", "max_intron"]:
            if key in rna and (
                not isinstance(rna[key], int)
                or isinstance(rna[key], bool)
                or rna[key] < (0 if key == "min_mapq" else 1)
            ):
                raise ValueError("Invalid RNA " + key)
        if rna.get("max_intron", 500000) < rna.get("min_intron", 20):
            raise ValueError("RNA max_intron must not be below min_intron")
        if "skip_duplicates" in rna and not isinstance(rna["skip_duplicates"], bool):
            raise ValueError("RNA skip_duplicates must be boolean")
    result["parameters"] = params
    for source, spec in result['predictions'].items():
        if 'residue_profiles' not in spec:
            continue
        try:
            import h5py
        except ImportError as exc:
            raise RuntimeError('Residue profiles require genesensus[molecular]') from exc
        with h5py.File(spec['residue_profiles'], 'r') as profiles:
            interpretation = str(profiles.attrs.get('interpretation', '')).lower()
        if 'binary-negative' in interpretation:
            if params['localizer_gate'] is None or 'scores' not in spec:
                raise ValueError(
                    'Residue profiles for '+source+' declare binary-negative scope; '
                    'provide binary scores and parameters.localizer_gate from the '
                    'independently validated score model. Positive/unclassified '
                    'profiles are excluded, not assigned zero risk.'
                )
    return result


def bam_index_paths(bam):
    bam = Path(bam)
    possible = [
        Path(str(bam) + ".bai"),
        bam.with_suffix(".bai"),
        Path(str(bam) + ".csi"),
        bam.with_suffix(".csi"),
    ]
    indexes = sorted({str(p.resolve()) for p in possible if p.is_file()})
    if not indexes:
        raise ValueError("RNA BAM index is missing")
    return indexes


def input_paths(config):
    out = {config["genome"]}
    for spec in config["predictions"].values():
        out.update(spec.values())
    for key in ["peptides", "rna", "proteins"]:
        for k, v in config.get(key, {}).items():
            if k in PATH_KEYS:
                out.add(v)
    if "bam" in config.get("rna", {}):
        out.update(bam_index_paths(config["rna"]["bam"]))
    return sorted(out)
