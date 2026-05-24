#!/usr/bin/env Rscript
# Generate R reference plots for SCORPIUS::draw_trajectory_plot.
suppressPackageStartupMessages({
  library(SCORPIUS)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Generate deterministic synthetic data
set.seed(42)
dataset <- generate_dataset(num_genes = 200, num_samples = 200, num_groups = 4)
space <- reduce_dimensionality(dataset$expression, ndim = 2)
groups <- dataset$sample_info$group_name
traj <- infer_trajectory(space)

# Save coordinates so Py side uses exactly the same numbers
write.csv(space, file.path(out_dir, "space.csv"), row.names = FALSE)
write.csv(data.frame(group = groups), file.path(out_dir, "groups.csv"), row.names = FALSE)
write.csv(traj$path, file.path(out_dir, "path.csv"), row.names = FALSE)
write.csv(data.frame(time = traj$time), file.path(out_dir, "time.csv"), row.names = FALSE)

# 1. Plain
p1 <- draw_trajectory_plot(space, progression_group = groups)
ggsave(file.path(out_dir, "R_traj_plain.png"), p1, width = 5, height = 5, dpi = 100)

# 2. With path
p2 <- draw_trajectory_plot(space, progression_group = groups, path = traj$path)
ggsave(file.path(out_dir, "R_traj_path.png"), p2, width = 5, height = 5, dpi = 100)

# 3. Numeric coloring
p3 <- draw_trajectory_plot(space, progression_group = traj$time)
ggsave(file.path(out_dir, "R_traj_numeric.png"), p3, width = 5, height = 5, dpi = 100)

cat("R plots:", file.path(out_dir, "R_traj_*.png"), "\n")
