# Discovery — py-SCORPIUS

## 1. Is the target already ported?

`gh repo view omicverse/py-SCORPIUS` → **not found at port start** (this repo is the first port).

## 2. R dependencies + py-mirror reuse

R upstream `DESCRIPTION` lists: `ggplot2`, `MASS`, `pheatmap`, `princurve`, `randomForest`, `pbapply`, `mclust`, `dplyr`, `purrr`, `dynutils`.

| R dep | Already mirrored under `omicverse/`? | Reused as |
|---|---|---|
| ggplot2 | ✅ via `Bio-Babel/ggplot2-python` | `ggplot2-python` (optional, for plotting) |
| pheatmap | ✅ via `Bio-Babel/pheatmap-python` | `pheatmap-python` (optional, for heatmaps) |
| mclust | ✅ `omicverse/py-mclustR` v0.2.0 | direct dep `pymclustR>=0.2.0` for `extract_modules` |
| princurve | ❌ no Python mirror | implemented inline (`pyscorpius.trajectory._principal_curve` — LOWESS-based Hastie-Stuetzle) |
| randomForest | n/a | sklearn `RandomForestRegressor` |
| MASS, pbapply, dynutils | n/a (utility) | numpy / scipy / tqdm |

## 3. Decision

**Proceed with full port** — algorithm class is ordinal (linear pseudotime), no upstream Python parity exists. Reuses `py-mclustR` for the Mclust module-extraction step.
