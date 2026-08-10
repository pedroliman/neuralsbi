#!/usr/bin/env Rscript
# Train neuralsbi's NRE on the shared simulations and sample each observation set.
# Usage: Rscript 11_run_neuralsbi_nre.R --classifier resnet --n_samples 5000
suppressMessages(library(neuralsbi))

args <- as.list(commandArgs(trailingOnly = TRUE))
opt <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i)) default else args[[i + 1L]]
}
classifier <- opt("--classifier", "resnet")
n_samples <- as.integer(opt("--n_samples", "5000"))
prior_var <- as.numeric(opt("--prior_var", "0.1"))
noise_var <- as.numeric(opt("--noise_var", "0.1"))
seed <- as.integer(opt("--seed", "42"))

data_dir <- file.path("data", "gaussian_linear_iid_nre")
meta <- as.integer(readLines(file.path(data_dir, "meta.txt")))
dim <- meta[1]
grid_n <- meta[length(meta)]
n_obs <- meta[2:(length(meta) - 1)]

theta <- as.matrix(read.csv(file.path(data_dir, "theta.csv"), header = FALSE))
x <- as.matrix(read.csv(file.path(data_dir, "x.csv"), header = FALSE))
task <- task_gaussian_linear(dim = dim, prior_var = prior_var, noise_var = noise_var)

fit <- nre(task$prior, theta = theta, x = x,
           classifier = classifier, seed = seed, verbose = TRUE)

out_dir <- file.path("results", "gaussian_linear_iid_nre")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (k in n_obs) {
  x_obs <- as.matrix(read.csv(sprintf("%s/x_obs_n%d.csv", data_dir, k),
                              header = FALSE))
  post <- posterior(fit, x_obs, seed = seed)
  draws <- sample(post, n_samples)
  cat(sprintf("n = %d: max rhat %.3f, min bulk ESS %.0f\n", k,
              max(attr(draws, "diagnostics")$rhat),
              min(attr(draws, "diagnostics")$ess_bulk)))
  path <- file.path(out_dir, sprintf("neuralsbi_%s_n%d.csv", classifier, k))
  write.table(unclass(draws), path, sep = ",", row.names = FALSE,
              col.names = FALSE)
  cat("wrote", path, "\n")
}

# Secondary metric (#146): the learned ratio itself, on the grid
# 09_generate_data_nre.R wrote around the true parameter, conditioned on the
# first observation.
grid_path <- file.path(data_dir, "grid_theta.csv")
x_obs1_path <- file.path(data_dir, "x_obs_n1.csv")
if (grid_n > 0 && file.exists(grid_path) && file.exists(x_obs1_path)) {
  grid_theta <- as.matrix(read.csv(grid_path, header = FALSE))
  x_obs_1 <- as.matrix(read.csv(x_obs1_path, header = FALSE))
  grid_ratio <- log_ratio(fit, grid_theta, x_obs_1)
  path <- file.path(out_dir,
                    sprintf("grid_log_ratio_neuralsbi_%s.csv", classifier))
  write.table(grid_ratio, path, sep = ",", row.names = FALSE, col.names = FALSE)
  cat("wrote", path, "\n")
}
