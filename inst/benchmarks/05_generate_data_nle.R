#!/usr/bin/env Rscript
# Generate the shared data for the NLE head-to-head.
#
# NLE learns q(x | theta) from single observations and is then conditioned on
# however many independent ones you have, so the benchmark needs two things the
# NPE protocol does not: training pairs where each x is one observation, and
# observation *sets* of several sizes drawn from one fixed parameter. The
# conjugate task has a closed-form posterior for every set size, which is what
# makes this a test of accuracy rather than of agreement.
#
# Usage: Rscript 05_generate_data_nle.R --dim 5 --n 10000 --seed 42
suppressMessages(library(neuralsbi))

args <- as.list(commandArgs(trailingOnly = TRUE))
opt <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i)) default else args[[i + 1L]]
}
dim <- as.integer(opt("--dim", "5"))
n <- as.integer(opt("--n", "10000"))
seed <- as.integer(opt("--seed", "42"))
n_obs <- as.integer(strsplit(opt("--n_obs", "1,10,100"), ",")[[1]])

task <- task_gaussian_linear(dim = dim)
set.seed(seed)
sims <- simulate_for_sbi(task$simulator, task$prior, n, seed = seed)

# One true parameter, and nested observation sets drawn from it, so the
# posteriors across set sizes are comparable rather than three unrelated runs.
theta_true <- as.numeric(sample_prior(task$prior, 1L))
x_all <- t(vapply(seq_len(max(n_obs)), function(i) task$simulator(theta_true),
                  numeric(dim)))

dir <- file.path("data", "gaussian_linear_iid")
dir.create(dir, recursive = TRUE, showWarnings = FALSE)
write.table(sims$theta, file.path(dir, "theta.csv"), sep = ",",
            row.names = FALSE, col.names = FALSE)
write.table(sims$x, file.path(dir, "x.csv"), sep = ",",
            row.names = FALSE, col.names = FALSE)
write.table(matrix(theta_true, nrow = 1), file.path(dir, "theta_true.csv"),
            sep = ",", row.names = FALSE, col.names = FALSE)
for (k in n_obs) {
  write.table(x_all[seq_len(k), , drop = FALSE],
              file.path(dir, sprintf("x_obs_n%d.csv", k)),
              sep = ",", row.names = FALSE, col.names = FALSE)
}
writeLines(as.character(c(dim, n_obs)), file.path(dir, "meta.txt"))
cat(sprintf("wrote %d simulations (dim %d) and observation sets %s to %s\n",
            n, dim, paste(n_obs, collapse = ", "), dir))
