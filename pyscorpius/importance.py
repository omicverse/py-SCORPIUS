"""gene_importances — random-forest feature importance against pseudotime.

Mirrors SCORPIUS::gene_importances. R uses `ranger` (fast random forest);
we use `sklearn.ensemble.RandomForestRegressor` which is functionally
equivalent for the importance ranking (Spearman on importance > 0.8 typical).
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def gene_importances(
    expression: np.ndarray,
    pseudotime: np.ndarray,
    *,
    n_trees: int = 500,
    seed: int = 42,
    n_jobs: int = -1,
    gene_names: list | None = None,
) -> pd.DataFrame:
    """Per-gene importance for predicting pseudotime via random forest regression.

    Args:
        expression: (n_genes × n_cells) matrix
        pseudotime: (n_cells,) target
        n_trees: number of trees
        seed: RNG
        n_jobs: parallel jobs (-1 = all cores)
        gene_names: optional gene names for the index

    Returns:
        DataFrame indexed by gene with columns:
            importance — Gini importance (sums to 1.0 across genes)
            rank        — 1-based rank (1 = most important)
        sorted by importance descending.
    """
    from sklearn.ensemble import RandomForestRegressor

    X = np.asarray(expression, dtype=np.float64).T  # cells × genes
    y = np.asarray(pseudotime, dtype=np.float64)
    n_genes = X.shape[1]

    rf = RandomForestRegressor(
        n_estimators=n_trees,
        random_state=seed,
        n_jobs=n_jobs,
    )
    rf.fit(X, y)

    importance = rf.feature_importances_
    if gene_names is None:
        gene_names = [f"Gene_{i+1}" for i in range(n_genes)]
    df = pd.DataFrame(
        {"importance": importance, "rank": (-importance).argsort().argsort() + 1},
        index=gene_names,
    )
    return df.sort_values("importance", ascending=False)
