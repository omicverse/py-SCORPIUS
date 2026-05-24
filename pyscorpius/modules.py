"""extract_modules — gene-module clustering via Mclust on expression-vs-pseudotime.

Mirrors SCORPIUS::extract_modules. R uses mclust on smoothed expression
patterns; we use py-mclustR (≥ 0.2.0) for the same.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d


def _smooth_along_time(expression: np.ndarray, pseudotime: np.ndarray, n_bins: int = 50,
                       sigma: float = 1.5) -> np.ndarray:
    """Smooth expression along pseudotime: bin → mean → Gaussian smooth.

    Returns (n_genes × n_bins) matrix of smoothed expression patterns.
    """
    n_genes, n_cells = expression.shape
    bins = np.linspace(0, 1, n_bins + 1)
    order = np.argsort(pseudotime)
    pt_sorted = pseudotime[order]
    expr_sorted = expression[:, order]
    binned = np.zeros((n_genes, n_bins))
    for b in range(n_bins):
        mask = (pt_sorted >= bins[b]) & (pt_sorted < bins[b + 1])
        if mask.sum() == 0:
            # take nearest non-empty bin
            binned[:, b] = expr_sorted[:, max(0, b - 1)] if b > 0 else expr_sorted[:, 0]
        else:
            binned[:, b] = expr_sorted[:, mask].mean(axis=1)
    # Gaussian-smooth each gene's binned trajectory
    smoothed = gaussian_filter1d(binned, sigma=sigma, axis=1)
    return smoothed


def extract_modules(
    expression: np.ndarray,
    pseudotime: np.ndarray,
    *,
    n_modules_range: range = range(2, 9),
    n_bins: int = 50,
    sigma: float = 1.5,
    seed: int = 42,
) -> dict:
    """Cluster genes into modules by their smoothed expression along pseudotime.

    Args:
        expression: (n_genes × n_cells) matrix
        pseudotime: (n_cells,) in [0, 1]
        n_modules_range: candidate G values for Mclust BIC scan
        n_bins: number of bins along pseudotime for the smoothing step
        sigma: Gaussian smoothing scale (in bins)
        seed: RNG seed

    Returns:
        dict with keys:
            module: (n_genes,) integer module ID per gene
            n_modules: int — selected G
            smoothed: (n_genes × n_bins) — the smoothed patterns
            bic: BIC value of the selected fit
    """
    expression = np.asarray(expression, dtype=np.float64)
    pseudotime = np.asarray(pseudotime, dtype=np.float64)

    # Smooth each gene's expression along pseudotime
    smoothed = _smooth_along_time(expression, pseudotime, n_bins=n_bins, sigma=sigma)
    # Z-score per gene (so we cluster by *shape* not magnitude)
    smoothed_z = smoothed - smoothed.mean(axis=1, keepdims=True)
    sd = smoothed_z.std(axis=1, keepdims=True, ddof=1)
    sd[sd == 0] = 1.0
    smoothed_z = smoothed_z / sd

    # Cluster genes with Mclust; fall back to KMeans + BIC if Mclust fails to
    # converge (e.g., too few genes relative to bin dimensionality)
    Gs = [g for g in n_modules_range if g > 1]
    try:
        from mclust_py import Mclust
        res = Mclust(smoothed_z, G=Gs, model_names=["VVV"])
        if res is None or not hasattr(res, 'classification'):
            raise RuntimeError("Mclust returned no fit")
        module_ids = res.classification.astype(int)
        if module_ids.min() == 0:
            module_ids = module_ids + 1
        return {
            "module": module_ids,
            "n_modules": int(res.G),
            "smoothed": smoothed,
            "bic": float(res.bic),
            "method": "mclust",
        }
    except (ImportError, RuntimeError, Exception):
        # Fallback: KMeans + BIC-style selection
        from sklearn.cluster import KMeans
        from sklearn.metrics import silhouette_score
        best_g, best_score = Gs[0], -np.inf
        best_km = None
        for g in Gs:
            try:
                km = KMeans(n_clusters=g, random_state=seed, n_init=10).fit(smoothed_z)
                if g < smoothed_z.shape[0]:
                    score = silhouette_score(smoothed_z, km.labels_)
                else:
                    score = -np.inf
                if score > best_score:
                    best_score = score; best_g = g; best_km = km
            except Exception:
                continue
        if best_km is None:
            best_km = KMeans(n_clusters=Gs[0], random_state=seed, n_init=10).fit(smoothed_z)
            best_g = Gs[0]
        module_ids = best_km.labels_ + 1
        return {
            "module": module_ids,
            "n_modules": int(best_g),
            "smoothed": smoothed,
            "bic": float("nan"),
            "method": "kmeans-fallback",
        }
