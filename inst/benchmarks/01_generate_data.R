#!/usr/bin/env Rscript
# Generate shared simulations for a benchmark task.
# Usage: Rscript 01_generate_data.R --task gaussian_linear --n 10000 --seed 42
suppressMessages(library(neuralsbi))

args <- as.list(commandArgs(trailingOnly = TRUE))
opt <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i)) default else args[[i + 1L]]
}
task_name <- opt("--task", "gaussian_linear")
n <- as.integer(opt("--n", "10000"))
n_obs <- as.integer(opt("--n_obs", "5"))
seed <- as.integer(opt("--seed", "42"))

task <- switch(task_name,
  gaussian_linear = task_gaussian_linear(),
  two_moons = task_two_moons(),
  slcp = task_slcp(),
  stop("unknown task: ", task_name)
)

set.seed(seed)
sims <- simulate_for_sbi(task$simulator, task$prior, n, seed = seed)
# observations: simulate from fresh prior draws (kept for reference too)
theta_obs <- sample_prior(task$prior, n_obs)
x_obs <- task$simulator(theta_obs)

dir <- file.path("data", task_name)
dir.create(dir, recursive = TRUE, showWarnings = FALSE)
write.table(sims$theta, file.path(dir, "theta.csv"), sep = ",",
            row.names = FALSE, col.names = FALSE)
write.table(sims$x, file.path(dir, "x.csv"), sep = ",",
            row.names = FALSE, col.names = FALSE)
write.table(theta_obs, file.path(dir, "theta_obs.csv"), sep = ",",
            row.names = FALSE, col.names = FALSE)
write.table(x_obs, file.path(dir, "x_obs.csv"), sep = ",",
            row.names = FALSE, col.names = FALSE)
cat(sprintf("wrote %d simulations + %d observations to %s\n", n, n_obs, dir))
