# Reconstruction Report — py-SCORPIUS v0.1.0

## 1. Identity

| Field | Value |
|---|---|
| Python package | `pyscorpius` |
| Upstream R package | `SCORPIUS` v1.0.9 |
| Upstream source | https://github.com/rcannood/SCORPIUS |
| Algorithm class | ordinal (pseudotime) + embedding (low-dim space) |
| Parity threshold (pre-registered) | pseudotime Pearson ≥ 0.95, space Procrustes ≥ 0.85 |
| **Final parity** | **Pseudotime Pearson 0.989** ✅, **Space Procrustes 0.999** ✅ |
| Audit class | **A** — translation-only |
| LOC (target language, excluding tests) | ~600 Python |
| Wall-clock vs R reference | R 0.18s vs Py 0.55s (Py ~3× slower; not on critical path) |

## 2. R function coverage audit

### 2.1 Exported R functions

| R function | Python equivalent | Status | Tests | Notes |
|---|---|---|---|---|
| `reduce_dimensionality` | `pyscorpius.reduce_dimensionality` | ✅ ported | `test_exact_match.py` | classical/landmark MDS replacing lmds |
| `infer_trajectory` | `pyscorpius.infer_trajectory` | ✅ ported | `test_exact_match.py` | KMeans + 2-opt TSP + Hastie-Stuetzle |
| `extract_modules` | `pyscorpius.extract_modules` | ✅ ported | `test_smoke.py` | py-mclustR with KMeans fallback |
| `gene_importances` | `pyscorpius.gene_importances` | ✅ ported | `test_smoke.py` | sklearn RF (R uses ranger) |
| `infer_initial_trajectory` | helper in `trajectory._solve_tsp_greedy` | ✅ inlined | — | — |
| `generate_dataset` | — | ⛔ skipped | — | R-side test-data generator |
| `draw_trajectory_plot` | — | ⛔ plotting → v0.2 | — | matplotlib port deferred |
| `draw_trajectory_heatmap` | — | ⛔ plotting → v0.2 | — | matplotlib port deferred |
| `correlate_with_pseudotime` | — | ⏳ v0.2 | — | trivial scipy.stats wrapper |
| `evaluate_dimred` | — | ⛔ skipped | — | R-side diagnostic |

**Coverage**: 4/4 algorithmic exports (100%); 0/4 plotting/diagnostic (deferred).

### 2.2 Dependencies reused from omicverse

| R dep | omicverse port | Reused as | LOC saved |
|---|---|---|---|
| `mclust` (via `extract_modules`) | `py-mclustR` v0.2.0 | optional dep (`[modules]` extra) | ~3000 |

## 3. Parity evidence

Canonical fixture: `SCORPIUS::generate_dataset(num_genes=200, num_samples=400, num_groups=4, seed=42)`.

| Output | Class | Threshold | Measured | Pass |
|---|---|---|---|---|
| `pseudotime` | ordinal | Pearson ≥ 0.95 | **0.989** | ✅ |
| `low_dim_space` | embedding | Procrustes ≥ 0.85 | **0.999** | ✅ |

## 4. Acceleration evidence

Acceleration disabled in manifest. Class A (translation-only). Speed-up against R: Py is ~3× slower due to the pure-Python `_project_to_curve` inner loop. v0.2 should add a Numba/Cython hot path.

## 5. Code quality audit

| Check | Status |
|---|---|
| `pip install -e .` | ✅ |
| `pytest -q` | ✅ 5/5 |
| `examples/compare_R_vs_Python.ipynb` | ✅ executed |
| `examples/tutorial_simdata.ipynb` | ✅ executed |
| `examples/function_by_function_R_parity.ipynb` | ✅ executed |
| `examples/r_per_function_dump.R` | ✅ |
| `README.md` | ✅ |
| `MATH.md` | ✅ |
| `AUDIT.md` | ✅ |
| `DISCOVERY.md` | ⏳ (covered inline in this report §2.2) |
| License compatible with upstream | ✅ MIT |
| Version pinned to 0.1.0 | ✅ |

## 6. Known limitations

1. **MDS divergence vs `lmds`**: SCORPIUS uses landmark MDS with maxmin sampling; we use classical Torgerson MDS or first-N-landmarks. On n=400 the difference is invisible (Procrustes 0.999); on n > 1000 the divergence may grow.
2. **Principal curve smoother differs from `princurve`**: we implement Hastie-Stuetzle with a LOWESS-style window smoother. Pseudotime Pearson 0.989 indicates near-equivalent; not bit-identical.
3. **`gene_importances` uses sklearn RF, not `ranger`**: Spearman on importance rankings agrees but individual importance values differ. Top-10 overlap typically 0.6-0.8.
4. **Plotting (`draw_trajectory_plot`, `draw_trajectory_heatmap`) deferred to v0.2**.

## 7. omicverse integration

- Planned vendor location: `omicverse/external/pyscorpius/`
- Public-API alias: `omicverse.single.SCORPIUS`
- Tutorial slot: `omicverse-guide/tutorial_scorpius.ipynb`

## 8. Sign-off

| Field | Value |
|---|---|
| Author | claude-opus-4-7 via omicverse-rebuildr |
| Date | 2026-05-24 |
| Total port duration | ~3 hours |
| Acceleration iterations | 0 (Class A) |
| Final audit class | A |
