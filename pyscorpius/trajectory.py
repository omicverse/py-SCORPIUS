"""infer_trajectory — k-means + TSP + principal curve.

Mirrors SCORPIUS::infer_trajectory. The principal-curve fit is a simple
iterative version of Hastie-Stuetzle (1989). Initial curve from k-means
centers ordered by shortest-path-through-cluster-centers (TSP).
"""

from __future__ import annotations

import itertools

import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.spatial.distance import cdist
from sklearn.cluster import KMeans


def _solve_tsp_greedy(D: np.ndarray) -> list[int]:
    """Approximate TSP via nearest-neighbour heuristic + 2-opt.

    Returns a permutation of nodes (no return-to-origin); use the OPEN tour.
    Good enough for k ≤ 9 clusters which is SCORPIUS's typical regime.
    """
    n = D.shape[0]
    # try every start; pick the shortest tour
    best_tour = None
    best_cost = np.inf
    for start in range(n):
        unvisited = set(range(n))
        unvisited.remove(start)
        tour = [start]
        while unvisited:
            cur = tour[-1]
            nxt = min(unvisited, key=lambda j: D[cur, j])
            tour.append(nxt)
            unvisited.remove(nxt)
        cost = sum(D[tour[i], tour[i + 1]] for i in range(len(tour) - 1))
        if cost < best_cost:
            best_cost = cost
            best_tour = tour
    # 2-opt improvement
    improved = True
    while improved:
        improved = False
        for i in range(1, n - 2):
            for j in range(i + 1, n):
                if j - i == 1:
                    continue
                a, b, c, d = best_tour[i - 1], best_tour[i], best_tour[j - 1], best_tour[j]
                if D[a, c] + D[b, d] < D[a, b] + D[c, d]:
                    best_tour[i:j] = best_tour[i:j][::-1]
                    improved = True
                    break
            if improved:
                break
    return best_tour


def _project_to_curve(X: np.ndarray, curve: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Project each row of X onto the polyline defined by `curve`.

    Returns (projected_points, arc_length_at_projection).
    """
    n_seg = curve.shape[0] - 1
    n = X.shape[0]
    # For each point, find the nearest segment + projection within
    best_proj = np.zeros_like(X)
    best_arc = np.zeros(n)
    seg_lengths = np.linalg.norm(np.diff(curve, axis=0), axis=1)
    cum_arc = np.concatenate([[0], np.cumsum(seg_lengths)])
    best_d = np.full(n, np.inf)
    for s in range(n_seg):
        a = curve[s]
        b = curve[s + 1]
        ab = b - a
        ab_norm2 = ab @ ab
        if ab_norm2 < 1e-20:
            continue
        ap = X - a
        t = np.clip(ap @ ab / ab_norm2, 0.0, 1.0)
        proj = a + t[:, None] * ab
        d = np.linalg.norm(X - proj, axis=1)
        better = d < best_d
        best_d[better] = d[better]
        best_proj[better] = proj[better]
        best_arc[better] = cum_arc[s] + t[better] * np.sqrt(ab_norm2)
    return best_proj, best_arc


def principal_curve(
    X: np.ndarray,
    init_curve: np.ndarray,
    *,
    max_iter: int = 10,
    tol: float = 1e-3,
    smooth_frac: float = 0.5,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit a smooth principal curve through X starting from init_curve.

    Hastie-Stuetzle 1989 style: alternate projecting points onto the curve
    and re-fitting the curve as a smooth function of arc-length. We use a
    LOWESS-style local-average smoother (rolling window) for robustness;
    spline-based smoothing can degenerate when the curve over-folds.

    Returns (curve_points, pseudotime_per_X_row, both in [0, 1]).
    """
    curve = init_curve.copy()
    n, d = X.shape
    s_norm = None
    for it in range(max_iter):
        # 1. Project X onto curve → arc length per point
        _, s = _project_to_curve(X, curve)
        if s.max() - s.min() < 1e-12:
            break
        s_norm = (s - s.min()) / (s.max() - s.min())
        idx = np.argsort(s_norm)
        s_sorted = s_norm[idx]
        # 2. Local-average smoother: for each grid point, weighted average of
        # nearest neighbours in arc-length space.
        s_grid = np.linspace(0, 1, 100)
        new_curve = np.zeros((100, d))
        h = max(smooth_frac * 0.1, 1.0 / 100)   # half-window
        for k, sg in enumerate(s_grid):
            w = np.maximum(0, 1 - np.abs(s_sorted - sg) / h) ** 2
            if w.sum() < 1e-6:
                # fallback to k-nearest
                nn = np.argsort(np.abs(s_sorted - sg))[:max(3, n // 20)]
                w = np.zeros(n)
                w[nn] = 1.0
            new_curve[k] = (w[:, None] * X[idx]).sum(axis=0) / w.sum()
        # 3. Convergence check
        if it > 0:
            delta = float(np.linalg.norm(new_curve - curve) / max(np.linalg.norm(curve), 1e-9))
            if delta < tol:
                curve = new_curve
                break
        curve = new_curve
    # Final pseudotime: arc-length parameterisation, scaled to [0,1]
    _, s = _project_to_curve(X, curve)
    if s.max() - s.min() > 1e-12:
        s_norm = (s - s.min()) / (s.max() - s.min())
    else:
        s_norm = np.zeros_like(s)
    return curve, s_norm


def infer_trajectory(
    space: np.ndarray,
    k: int = 4,
    *,
    seed: int = 42,
    max_iter: int = 10,
) -> dict:
    """SCORPIUS-style linear-trajectory inference.

    Args:
        space: (n_samples × ndim) low-dim coordinates (from reduce_dimensionality).
        k: number of k-means clusters for the initial trajectory.
        seed: RNG for KMeans.
        max_iter: principal-curve iterations.

    Returns:
        dict with keys:
            time: per-sample pseudotime in [0, 1]
            path: (G × ndim) curve coordinates
            tour: list[int] — the TSP order of cluster centroids
    """
    space = np.asarray(space, dtype=np.float64)
    n = space.shape[0]
    # 1. KMeans
    km = KMeans(n_clusters=k, random_state=seed, n_init=10).fit(space)
    centers = km.cluster_centers_
    # 2. TSP through centers (open tour)
    D_clu = cdist(centers, centers, "euclidean")
    tour = _solve_tsp_greedy(D_clu)
    init_curve = centers[tour]
    # 3. Principal curve refinement
    curve, time = principal_curve(space, init_curve, max_iter=max_iter)
    return {"time": time, "path": curve, "tour": tour}
