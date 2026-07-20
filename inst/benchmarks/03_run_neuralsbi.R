#!/usr/bin/env Rscript
# Train neuralsbi's NPE on the shared simulations and save posterior samples.
# Usage: Rscript 03_run_neuralsbi.R --task gaussian_linear --estimator maf
suppressMessages(library(neuralsbi))

args <- as.list(commandArgs(trailingOnly = TRUE))
opt <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i)) default else args[[i + 1L]]
}
task_name <- opt("--task", "gaussian_linear")
estimator <- opt("--estimator", "maf")
n_samples <- as.integer(opt("--n_samples", "10000"))
seed <- as.integer(opt("--seed", "42"))

task <- switch(task_name,
  gaussian_linear = task_gaussian_linear(),
  two_moons = task_two_moons(),
  slcp = task_slcp(),
  stop("unknown task: ", task_name)
)

data_dir <- file.path("data", task_name)
theta <- as.matrix(read.csv(file.path(data_dir, "theta.csv"), header = FALSE))
x <- as.matrix(read.csv(file.path(data_dir, "x.csv"), header = FALSE))
x_obs <- as.matrix(read.csv(file.path(data_dir, "x_obs.csv"), header = FALSE))

fit <- npe(task$prior, theta = theta, x = x,
           density_estimator = estimator, seed = seed, verbose = TRUE)

out_dir <- file.path("results", task_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (i in seq_len(nrow(x_obs))) {
  post <- posterior(fit, x_obs = x_obs[i, ])
  draws <- sample(post, n_samples)
  path <- file.path(out_dir, sprintf("neuralsbi_%s_obs%d.csv", estimator, i))
  write.table(unclass(draws), path, sep = ",", row.names = FALSE,
              col.names = FALSE)
  cat("wrote", path, "\n")
}
