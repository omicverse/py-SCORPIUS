#!/usr/bin/env Rscript
# Per-function dump for SCORPIUS Notebook 3.
suppressMessages({
  library(SCORPIUS); library(jsonlite)
})
OUT <- "examples/_r_outputs"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

set.seed(42)
ds <- SCORPIUS::generate_dataset(num_genes = 200, num_samples = 400, num_groups = 4)
expression <- ds$expression
group <- ds$sample_info$group_name

set.seed(42)
space <- SCORPIUS::reduce_dimensionality(expression, dist = "spearman", ndim = 2)

set.seed(42)
traj <- SCORPIUS::infer_trajectory(space)

# Reduce dimensionality
write_json(list(
  args = list(dist = "spearman", ndim = 2, num_landmarks = 1000),
  space = as.matrix(space),
  shape = dim(space)
), file.path(OUT, "reduce_dimensionality.json"),
auto_unbox = TRUE, digits = NA, matrix = "rowmajor", na = "null", pretty = FALSE)

# Infer trajectory
write_json(list(
  args = list(k = 4),
  time = as.numeric(traj$time),
  path = as.matrix(traj$path),
  path_shape = dim(traj$path)
), file.path(OUT, "infer_trajectory.json"),
auto_unbox = TRUE, digits = NA, matrix = "rowmajor", na = "null", pretty = FALSE)

# Gene importances
set.seed(42)
imp <- SCORPIUS::gene_importances(expression, traj$time, num_permutations = 0)
write_json(list(
  args = list(num_permutations = 0),
  gene_names = rownames(imp),
  importance = as.numeric(imp$importance),
  top10 = head(rownames(imp), 10)
), file.path(OUT, "gene_importances.json"),
auto_unbox = TRUE, digits = NA, na = "null", pretty = FALSE)

cat("[r-dump] all done.\n")
