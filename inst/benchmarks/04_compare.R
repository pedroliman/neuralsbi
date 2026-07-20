#!/usr/bin/env Rscript
# Compare neuralsbi vs Python sbi posterior samples with C2ST and moments.
# Usage: Rscript 04_compare.R --task gaussian_linear --estimator maf
suppressMessages(library(neuralsbi))

args <- as.list(commandArgs(trailingOnly = TRUE))
opt <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i)) default else args[[i + 1L]]
}
task_name <- opt("--task", "gaussian_linear")
estimator <- opt("--estimator", "maf")

task <- switch(task_name,
  gaussian_linear = task_gaussian_linear(),
  two_moons = task_two_moons(),
  slcp = task_slcp(),
  stop("unknown task: ", task_name)
)
res_dir <- file.path("results", task_name)
x_obs <- as.matrix(read.csv(file.path("data", task_name, "x_obs.csv"),
                            header = FALSE))

rows <- list()
for (i in seq_len(nrow(x_obs))) {
  f_ours <- file.path(res_dir, sprintf("neuralsbi_%s_obs%d.csv", estimator, i))
  f_sbi <- file.path(res_dir, sprintf("sbi_%s_obs%d.csv", estimator, i))
  if (!file.exists(f_ours) || !file.exists(f_sbi)) next
  ours <- as.matrix(read.csv(f_ours, header = FALSE))
  theirs <- as.matrix(read.csv(f_sbi, header = FALSE))
  acc <- c2st(ours, theirs, seed = i)$accuracy
  mean_diff <- max(abs(colMeans(ours) - colMeans(theirs)))
  sd_diff <- max(abs(apply(ours, 2, sd) - apply(theirs, 2, sd)))
  row <- data.frame(obs = i, c2st_vs_sbi = acc,
                    max_mean_diff = mean_diff, max_sd_diff = sd_diff)
  if (!is.null(task$reference_posterior)) {
    ref <- task$reference_posterior(x_obs[i, ], nrow(ours))
    row$c2st_ours_vs_ref <- c2st(ours, ref, seed = i)$accuracy
    row$c2st_sbi_vs_ref <- c2st(theirs, ref, seed = i)$accuracy
  }
  rows[[length(rows) + 1L]] <- row
}
tab <- do.call(rbind, rows)
print(tab)
out <- sprintf("comparison_%s_%s.csv", task_name, estimator)
write.csv(tab, file.path(res_dir, out), row.names = FALSE)
cat("wrote", file.path(res_dir, out), "\n")
cat(sprintf("PASS criterion: c2st_vs_sbi <= 0.60 -> %s\n",
            if (all(tab$c2st_vs_sbi <= 0.60)) "PASS" else "FAIL"))
