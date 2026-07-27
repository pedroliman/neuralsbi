#!/usr/bin/env Rscript
# Train neuralsbi's NLE on the shared simulations and sample each observation set.
# Usage: Rscript 07_run_neuralsbi_nle.R --estimator maf --n_samples 5000
suppressMessages(library(neuralsbi))

args <- as.list(commandArgs(trailingOnly = TRUE))
opt <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i)) default else args[[i + 1L]]
}
estimator <- opt("--estimator", "maf")
n_samples <- as.integer(opt("--n_samples", "5000"))
seed <- as.integer(opt("--seed", "42"))

data_dir <- file.path("data", "gaussian_linear_iid")
meta <- as.integer(readLines(file.path(data_dir, "meta.txt")))
dim <- meta[1]
n_obs <- meta[-1]

theta <- as.matrix(read.csv(file.path(data_dir, "theta.csv"), header = FALSE))
x <- as.matrix(read.csv(file.path(data_dir, "x.csv"), header = FALSE))
task <- task_gaussian_linear(dim = dim)

fit <- nle(task$prior, theta = theta, x = x,
           density_estimator = estimator, seed = seed, verbose = TRUE)

out_dir <- file.path("results", "gaussian_linear_iid")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (k in n_obs) {
  x_obs <- as.matrix(read.csv(sprintf("%s/x_obs_n%d.csv", data_dir, k),
                              header = FALSE))
  post <- posterior(fit, x_obs, seed = seed)
  draws <- sample(post, n_samples)
  cat(sprintf("n = %d: max rhat %.3f, min bulk ESS %.0f\n", k,
              max(attr(draws, "diagnostics")$rhat),
              min(attr(draws, "diagnostics")$ess_bulk)))
  path <- file.path(out_dir, sprintf("neuralsbi_%s_n%d.csv", estimator, k))
  write.table(unclass(draws), path, sep = ",", row.names = FALSE,
              col.names = FALSE)
  cat("wrote", path, "\n")
}
