"""Vectorized CPU aggregation and experimental, hardware-unvalidated CUDA scoring."""

import numpy as np


def weighted_scores(values, available, weights, backend="cpu"):
    v = np.asarray(values, dtype=np.float64)
    mask = np.asarray(available, dtype=bool)
    w = np.asarray(weights, dtype=np.float64)
    if v.ndim != 2 or v.shape != mask.shape or w.shape != (v.shape[1],):
        raise ValueError("Incompatible score matrix shapes")
    if (
        not np.all(np.isfinite(v))
        or np.any(v < 0)
        or np.any(v > 1)
        or not np.all(np.isfinite(w))
        or np.any(w < 0)
    ):
        raise ValueError("Scores must be in [0,1] and weights finite/nonnegative")
    if backend not in ("cpu", "auto", "cuda"):
        raise ValueError("Backend must be cpu, auto or cuda")
    # Sparse graph/search work is CPU-bound. Auto stays CPU until profiling proves otherwise.
    if backend in ("cpu", "auto"):
        active = mask * w
        den = active.sum(axis=1)
        return np.divide(
            (v * active).sum(axis=1), den, out=np.zeros(len(v)), where=den > 0
        )
    try:
        import torch
    except ImportError as exc:
        raise RuntimeError("CUDA scoring requires genesensus[gpu]") from exc
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but no CUDA device is available")
    with torch.inference_mode():
        x = torch.as_tensor(v, device="cuda", dtype=torch.float64)
        a = torch.as_tensor(mask, device="cuda") * torch.as_tensor(
            w, device="cuda", dtype=torch.float64
        )
        den = a.sum(1)
        result = torch.where(
            den > 0,
            (x * a).sum(1) / den.clamp_min(1e-300),
            torch.zeros(len(v), device="cuda", dtype=torch.float64),
        )
        return result.cpu().numpy()
