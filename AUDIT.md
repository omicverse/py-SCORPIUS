## R function coverage audit

### Coverage summary

| Category | Ported | Total | % |
|---|---|---|---|
| Exported R functions | 0 | 15 | 0.0% |
| Internal helpers (reachable) | 0 | 4 | 0.0% |

_Python package exposes 0 unique names._

### Exported R functions

| R function | Python equivalent | Status |
|---|---|---|
| `apply_minmax_scale` | `—` | ❌ MISSING |
| `apply_quantile_scale` | `—` | ❌ MISSING |
| `apply_uniform_scale` | `—` | ❌ MISSING |
| `draw_trajectory_heatmap` | `—` | ❌ MISSING |
| `draw_trajectory_plot` | `—` | ❌ MISSING |
| `extract_modules` | `—` | ❌ MISSING |
| `gene_importances` | `—` | ❌ MISSING |
| `generate_dataset` | `—` | ❌ MISSING |
| `infer_trajectory` | `—` | ❌ MISSING |
| `reduce_dimensionality` | `—` | ❌ MISSING |
| `reverse_trajectory` | `—` | ❌ MISSING |
| `scale_minmax` | `—` | ❌ MISSING |
| `scale_quantile` | `—` | ❌ MISSING |
| `scale_uniform` | `—` | ❌ MISSING |
| `ti_scorpius` | `—` | ❌ MISSING |

### Internal helpers reachable from exports

| R helper | File | Python equivalent | Status |
|---|---|---|---|
| `check_logical_vector` | `dummy_proofing.R` | `—` | 🔸 missing-or-inlined |
| `check_numeric_matrix` | `dummy_proofing.R` | `—` | 🔸 missing-or-inlined |
| `check_numeric_vector` | `dummy_proofing.R` | `—` | 🔸 missing-or-inlined |
| `infer_initial_trajectory` | `trajectory_inference.R` | `—` | 🔸 missing-or-inlined |
