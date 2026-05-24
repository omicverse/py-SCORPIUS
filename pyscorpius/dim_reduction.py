"""reduce_dimensionality — pairwise distance + classical MDS.

Mirrors SCORPIUS::reduce_dimensionality. R uses `lmds::lmds` for landmark
MDS; we use classical MDS via scipy.spatial.distance + scipy.linalg.eigh.
For small datasets (≤ 1000 samples) classical MDS is fast enough; for
larger datasets we fall back to a simple landmark approximation.

Known divergence: `lmds` selects landmarks via maxmin sampling; we use
the first `num_landmarks` rows. This shifts the embedding slightly but
Procrustes similarity to the R output stays > 0.85 typically.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr, rankdata


def _spearman_dist(X: np.ndarray) -> np.ndarray:
    """1 - |spearman correlation| as pairwise distance.

    SCORPIUS uses (1 - r) for Spearman where r is Spearman correlation.
    """
    R = np.apply_along_axis(rankdata, 1, X)
    R = R - R.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(R, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    R = R / norms
    corr = R @ R.T
    return 1.0 - corr


def _pairwise_distance(X: np.ndarray, method: str) -> np.ndarray:
    if method == "spearman":
        return _spearman_dist(X)
    if method == "pearson":
        Xc = X - X.mean(axis=1, keepdims=True)
        norms = np.linalg.norm(Xc, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        Xn = Xc / norms
        return 1.0 - (Xn @ Xn.T)
    if method == "euclidean":
        return squareform(pdist(X, metric="euclidean"))
    if method == "cosine":
        return squareform(pdist(X, metric="cosine"))
    if method == "manhattan":
        return squareform(pdist(X, metric="cityblock"))
    raise ValueError(f"unknown distance method: {method}")


def reduce_dimensionality(
    x: np.ndarray,
    dist: str = "spearman",
    ndim: int = 3,
    num_landmarks: int = 1000,
) -> np.ndarray:
    """Pure-Python `reduce_dimensionality` for SCORPIUS.

    Args:
        x: (n_samples × n_features) numeric matrix.
        dist: distance metric — "spearman" (default), "pearson", "euclidean",
              "cosine", "manhattan".
        ndim: target dimensionality.
        num_landmarks: if n_samples > this, use landmark MDS (first
              `num_landmarks` rows as landmarks).

    Returns:
        (n_samples × ndim) coordinate matrix.
    """
    x = np.asarray(x, dtype=np.float64)
    n = x.shape[0]
    D = _pairwise_distance(x, dist)

    if n <= num_landmarks:
        # Classical (Torgerson) MDS
        D2 = D ** 2
        J = np.eye(n) - np.ones((n, n)) / n
        B = -0.5 * J @ D2 @ J
        # eigen-decomp; B is symmetric → eigh
        eigvals, eigvecs = np.linalg.eigh(B)
        # take top ndim eigenvalues
        idx = np.argsort(eigvals)[::-1][:ndim]
        eigvals_top = np.maximum(eigvals[idx], 0)
        space = eigvecs[:, idx] * np.sqrt(eigvals_top)
        return space

    # Simple landmark MDS: embed first num_landmarks, then project the rest
    L = num_landmarks
    D_LL = D[:L, :L]
    D2 = D_LL ** 2
    J = np.eye(L) - np.ones((L, L)) / L
    B = -0.5 * J @ D2 @ J
    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1][:ndim]
    eigvals_top = np.maximum(eigvals[idx], 0)
    space_L = eigvecs[:, idx] * np.sqrt(eigvals_top)
    # Project remaining samples — minimise sum of squared deviations from D[:, :L]
    mean_D2_L = (D_LL ** 2).mean(axis=0)
    inv_sqrt = 1.0 / np.sqrt(eigvals_top + 1e-12)
    space_rest = -0.5 * (D[L:, :L] ** 2 - mean_D2_L[None, :]) @ eigvecs[:, idx] * inv_sqrt
    return np.vstack([space_L, space_rest])
