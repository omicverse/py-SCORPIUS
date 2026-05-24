# py-SCORPIUS

A **pure-Python port of [SCORPIUS](https://github.com/rcannood/SCORPIUS)** (Cannoodt et al., bioRxiv 2016) for linear-trajectory inference in single-cell RNA-seq.

- AnnData-compatible (cells × genes)
- **No `rpy2`** — pure NumPy / SciPy / scikit-learn
- Same function surface as the R workflow (`reduce_dimensionality` → `infer_trajectory` → `extract_modules` → `gene_importances`)
- **Pseudotime Pearson = 0.989** vs R SCORPIUS on canonical fixture
- **Low-dim space Procrustes = 0.999** vs R

## Install

```bash
pip install pyscorpius
# Optional — for extract_modules:
pip install pyscorpius[modules]
```

## Quick-start

```python
import numpy as np
from pyscorpius import reduce_dimensionality, infer_trajectory

# expression: (n_cells × n_genes) — or pass an AnnData and use .X.toarray()
space = reduce_dimensionality(expression, dist="spearman", ndim=3)
traj  = infer_trajectory(space, k=4)
pseudotime = traj["time"]    # in [0, 1] per cell
curve      = traj["path"]    # smooth curve through `space`
```

## Function map

| Python | R counterpart | Purpose |
|---|---|---|
| `reduce_dimensionality` | `reduce_dimensionality` | distance + MDS to low-dim space |
| `infer_trajectory` | `infer_trajectory` | kmeans + TSP + Hastie-Stuetzle principal curve |
| `principal_curve` | `princurve::principal_curve` | exposed for direct use |
| `extract_modules` | `extract_modules` | gene-module clustering via Mclust (needs `pyscorpius[modules]`) |
| `gene_importances` | `gene_importances` | random-forest feature importance against pseudotime |

## Reproducing R results

```bash
# Run R reference under your R conda env
Rscript tests/r_reference_driver.R data/fixture_simdata.rds data/reference_output.json

# Run Python candidate
python tests/_run_candidate.py data/fixture_simdata.rds data/candidate_output.json

# Compare (see compare_R_vs_Python.ipynb)
pytest tests/test_exact_match.py -v
```

Achieved on the SCORPIUS-bundled simulated dataset (400 cells × 200 genes):
- Pseudotime Pearson: **0.989** (threshold 0.95)
- Low-dim space Procrustes: **0.999** (threshold 0.85)

## Relationship to omicverse

Developed under the [omicverse-rebuildr](https://github.com/omicverse/omicverse-rebuildr) protocol. Reuses [py-mclustR](https://github.com/omicverse/py-mclustR) for the `extract_modules` Mclust step.

## Citation

> Cannoodt, R. et al. **SCORPIUS improves trajectory inference and identifies novel modules in dendritic cell development.** bioRxiv 079509 (2016).

## License

MIT — matches upstream SCORPIUS.
