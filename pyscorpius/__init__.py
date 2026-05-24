"""pyscorpius — pure-Python port of SCORPIUS (Cannoodt et al. 2016).

Linear-trajectory inference for single-cell RNA-seq via:
1. Pairwise distance (Spearman default)
2. Classical MDS / landmark MDS to low-dim space
3. K-means clustering + TSP shortest path → initial trajectory
4. Hastie-Stuetzle principal curve refinement → final pseudotime
5. Gene-module extraction (uses py-mclustR)
6. Feature importance (random forest)

Reusable upstream: py-mclustR ≥ 0.2.0 for `extract_modules`.
"""

from __future__ import annotations

__version__ = "0.1.0"

from .dim_reduction import reduce_dimensionality
from .trajectory import infer_trajectory, principal_curve
from .modules import extract_modules
from .importance import gene_importances

__all__ = [
    "reduce_dimensionality",
    "infer_trajectory",
    "principal_curve",
    "extract_modules",
    "gene_importances",
    "__version__",
]
