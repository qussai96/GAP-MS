"""Explicit optional protein scores/profiles; no correctness labels or fitting."""

from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
import csv, gzip
import numpy as np


def read_scores(path, id_column="protein_id", score_column="protein_likeness_score"):
    opener = gzip.open if str(path).endswith(".gz") else open
    result = {}
    with opener(path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if (
            not reader.fieldnames
            or id_column not in reader.fieldnames
            or score_column not in reader.fieldnames
        ):
            raise ValueError(
                f"Score table requires {id_column} and {score_column}: {path}"
            )
        for row in reader:
            ident = row[id_column]
            q = float(row[score_column])
            if not ident or not np.isfinite(q) or not 0 <= q <= 1:
                raise ValueError(f"Invalid molecular score: {ident}")
            if ident in result and result[ident] != q:
                raise ValueError(f"Conflicting duplicate molecular score: {ident}")
            result[ident] = q
    if not result:
        raise ValueError("Empty score table")
    return result


@dataclass
class RiskTrack:
    model: object
    probabilities: object
    q: object = None

    _prefix: object = field(default=None, repr=False, init=False)
    _mean: object = field(default=None, repr=False, init=False)

    def risk(self, intervals, gate=None):
        # A diagnostic profile intended only for binary-negative proteins is
        # unobserved outside that scope, not a measured zero-error profile.
        if gate is not None and (self.q is None or self.q >= gate):
            return None, False
        if tuple(intervals) == self.model.cds:
            if not len(self.probabilities):
                return None, False
            if self._mean is None:
                self._mean = float(np.mean(self.probabilities, dtype=np.float64))
            return self._mean, True
        total = 0.0
        count = 0
        prefix = 0
        m = self.model
        if self._prefix is None:
            self._prefix = np.concatenate(
                ([0.0], np.cumsum(self.probabilities, dtype=np.float64))
            )

        def integral(position):
            k, r = divmod(position, 3)
            return float(self._prefix[k]) * 3 + (
                float(self.probabilities[k]) * r if r else 0.0
            )

        for a, b in m.cds:
            for x, y in intervals:
                left, right = max(a, x), min(b, y)
                if left >= right:
                    continue
                if m.strand == "+":
                    start, end = (
                        prefix + left - a - m.phase,
                        prefix + right - a - m.phase,
                    )
                else:
                    start, end = (
                        m.length - prefix - (right - a) - m.phase,
                        m.length - prefix - (left - a) - m.phase,
                    )
                start = max(0, start)
                end = min(len(self.probabilities) * 3, end)
                if start < end:
                    total += integral(end) - integral(start)
                    count += end - start
            prefix += b - a
        return (total / count, True) if count else (None, False)


def attach_molecular(candidates, source_models, specs, genome):
    """Validate profile IDs and sequence correspondence, then attach measured data."""
    bykey = {c.model.key: c for c in candidates}
    report = {}
    for source, models in source_models.items():
        spec = specs.get(source, {})
        if not spec:
            continue
        qs = read_scores(spec["scores"]) if spec.get("scores") else {}
        profiles = {}
        sequences = {}
        if spec.get("residue_profiles"):
            try:
                import h5py
            except ImportError as exc:
                raise RuntimeError(
                    "Residue profiles require genesensus[molecular]"
                ) from exc
            with h5py.File(spec["residue_profiles"], "r") as h:
                idkey = next((k for k in ["protein_id", "protein_ids"] if k in h), None)
                valkey = next(
                    (
                        k
                        for k in ["localization_probability", "probabilities"]
                        if k in h
                    ),
                    None,
                )
                if idkey is None or valkey is None or "offsets" not in h:
                    raise ValueError(
                        "HDF5 requires protein_id, offsets and localization_probability"
                    )
                if (
                    h[idkey].ndim != 1
                    or h["offsets"].ndim != 1
                    or h["offsets"].dtype.kind not in "iu"
                    or h[valkey].ndim != 1
                ):
                    raise ValueError(
                        "Residue IDs, integer offsets and probabilities must be one-dimensional"
                    )
                ids = h[idkey].asstr()[:].tolist()
                offsets = np.asarray(h["offsets"][:], dtype=np.int64)
                values = np.asarray(h[valkey][:])
                if any(not ident for ident in ids):
                    raise ValueError("Empty residue profile identifier")
                if "sequence" in h and (
                    h["sequence"].ndim != 1 or len(h["sequence"]) != len(ids)
                ):
                    raise ValueError("Residue sequence count differs from identifiers")
                if len(ids) != len(set(ids)):
                    raise ValueError("Duplicate residue profile IDs")
                if (
                    len(offsets) != len(ids) + 1
                    or offsets[0] != 0
                    or offsets[-1] != len(values)
                    or np.any(np.diff(offsets) < 0)
                ):
                    raise ValueError("Invalid residue profile offsets")
                if (
                    not np.all(np.isfinite(values))
                    or np.any(values < 0)
                    or np.any(values > 1)
                ):
                    raise ValueError(
                        "Residue probabilities must be finite and in [0,1]"
                    )
                profiles = {
                    tid: values[offsets[i] : offsets[i + 1]]
                    for i, tid in enumerate(ids)
                }
                if "sequence" in h:
                    sequences = dict(zip(ids, h["sequence"].asstr()[:].tolist()))
        observed = {}
        nq = nr = 0
        for m in models:
            if (
                m.identifier in observed
                and observed[m.identifier] != m.key
                and (m.identifier in qs or m.identifier in profiles)
            ):
                raise ValueError(
                    f"Ambiguous molecular identifier reused across structures: {source}:{m.identifier}"
                )
            observed[m.identifier] = m.key
            c = bykey[m.key]
            if m.identifier in qs:
                c.gene_scores.append(qs[m.identifier])
                nq += 1
            if m.identifier in profiles:
                pp = profiles[m.identifier]
                _, _, translated = genome.assess(m)
                if m.identifier in sequences:
                    expected = sequences[m.identifier].rstrip("*").upper()
                    if translated.rstrip("*") != expected:
                        raise ValueError(
                            f"Profile/annotation translation mismatch: {source}:{m.identifier}"
                        )
                elif len(pp) not in (len(translated), len(translated.rstrip("*"))):
                    raise ValueError(
                        f"Profile/annotation length mismatch: {source}:{m.identifier}"
                    )
                if m.identifier in sequences and len(pp) != len(
                    sequences[m.identifier]
                ):
                    raise ValueError(
                        f"Profile/sequence length mismatch: {m.identifier}"
                    )
                c.risks.append(RiskTrack(m, pp, qs.get(m.identifier)))
                nr += 1
        report[source] = dict(
            models=len(models),
            measured_scores=nq,
            measured_profiles=nr,
            unmatched_model_scores=len(models) - nq if qs else None,
            unmatched_model_profiles=len(models) - nr if profiles else None,
        )
        if qs and nq == 0:
            raise ValueError(f"No molecular score identifiers match source {source}")
        if profiles and nr == 0:
            raise ValueError(f"No residue identifiers match source {source}")
    return report


def local_risk(candidate, regions, gate=None):
    measured = [track.risk(regions, gate) for track in candidate.risks]
    values = [r for r, known in measured if known]
    return (float(np.mean(values)), True) if values else (None, False)
