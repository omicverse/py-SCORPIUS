#!/usr/bin/env Rscript
# SCORPIUS R reference: generate dataset + reduce_dimensionality + infer_trajectory
suppressMessages({
  library(SCORPIUS)
  library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
fixture_path <- args[1]; output_path <- args[2]

set.seed(42)
ds <- SCORPIUS::generate_dataset(num_genes = 200, num_samples = 400, num_groups = 4)
expression <- ds$expression
group <- ds$sample_info$group_name

# Save fixture for Python side
saveRDS(list(expression = expression, group = as.character(group)), fixture_path)
write.csv(expression, sub("\\.rds$", "_expression.csv", fixture_path))
write.csv(data.frame(group = group), sub("\\.rds$", "_group.csv", fixture_path),
          row.names = FALSE)

set.seed(42)
t0 <- proc.time()[["elapsed"]]
space <- SCORPIUS::reduce_dimensionality(expression, dist = "spearman", ndim = 2)
t_rd <- proc.time()[["elapsed"]] - t0
cat("[ref] reduce_dimensionality:", t_rd, "s; space dims:", dim(space), "\n")

set.seed(42)
t0 <- proc.time()[["elapsed"]]
traj <- SCORPIUS::infer_trajectory(space)
t_it <- proc.time()[["elapsed"]] - t0
cat("[ref] infer_trajectory:", t_it, "s; pseudotime length:", length(traj$time), "\n")

out <- list(
  cell_names = rownames(expression),
  group = as.character(group),
  space = as.matrix(space),
  pseudotime = as.numeric(traj$time),
  path = as.matrix(traj$path),
  timings = list(reduce_dimensionality = t_rd, infer_trajectory = t_it)
)
write_json(out, output_path, auto_unbox = TRUE, digits = NA, matrix = "rowmajor",
           na = "null", pretty = FALSE)
cat("[ref] wrote", output_path, "\n")
