#!/usr/bin/env Rscript
# Generate the shared data for the NRE head-to-head.
#
# nre() learns r(theta, x) = p(x | theta) / p(x) from single observations and
# conditions on however many independent ones it is given, so this needs the
# same observation-set layout 05_generate_data_nle.R uses for NLE. It adds one
# thing NLE does not need: a fixed grid of theta values, conditioned on the
# first observation, so 12_compare_nre.R can score the learned ratio itself
# against the analytic one -- not just the posteriors it implies. The grid is
# only worth building in low dimension (see #146's own 2-dimensional run), so
# it is skipped for dim > 2 and the metric that depends on it is skipped too.
#
# Usage: Rscript 09_generate_data_nre.R --dim 2 --n 8000 --seed 42
suppressMessages(library(neuralsbi))

args <- as.list(commandArgs(trailingOnly = TRUE))
opt <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i)) default else args[[i + 1L]]
}
dim <- as.integer(opt("--dim", "2"))
n <- as.integer(opt("--n", "8000"))
seed <- as.integer(opt("--seed", "42"))
n_obs <- as.integer(strsplit(opt("--n_obs", "1,20"), ",")[[1]])
prior_var <- as.numeric(opt("--prior_var", "0.1"))
noise_var <- as.numeric(opt("--noise_var", "0.1"))
grid_n <- as.integer(opt("--grid_n", "11"))
grid_width <- as.numeric(opt("--grid_width", "4"))

task <- task_gaussian_linear(dim = dim, prior_var = prior_var, noise_var = noise_var)
set.seed(seed)
sims <- simulate_for_sbi(task$simulator, task$prior, n, seed = seed)

# One true parameter, and nested observation sets drawn from it, exactly as
# 05_generate_data_nle.R does, so the posteriors across set sizes are
# comparable rather than three unrelated runs.
theta_true <- as.numeric(sample_prior(task$prior, 1L))
x_all <- t(vapply(seq_len(max(n_obs)), function(i) task$simulator(theta_true),
                  numeric(dim)))

dir <- file.path("data", "gaussian_linear_iid_nre")
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

# The grid conditions on the first observation, so it needs 1 in n_obs.
has_grid <- dim <= 2L && 1L %in% n_obs
if (has_grid) {
  # Half-width in units of the single-observation conjugate posterior sd, so
  # the grid covers where the ratio actually varies rather than an arbitrary
  # box (see task_gaussian_linear()'s reference() for the same formula).
  post_var <- 1 / (1 / prior_var + 1 / noise_var)
  half_width <- grid_width * sqrt(post_var)
  axes <- lapply(seq_len(dim), function(j)
    seq(theta_true[j] - half_width, theta_true[j] + half_width,
        length.out = grid_n))
  grid_theta <- as.matrix(do.call(expand.grid, axes))
  dimnames(grid_theta) <- NULL
  write.table(grid_theta, file.path(dir, "grid_theta.csv"), sep = ",",
              row.names = FALSE, col.names = FALSE)
} else if (dim > 2L) {
  cat(sprintf("dim = %d > 2: skipping the analytic-log-ratio grid.\n", dim))
}

# All-integer, one value per line, like 05_generate_data_nle.R's meta.txt:
# dim, then the n_obs values, then grid_n (0 if no grid was written).
writeLines(as.character(c(dim, n_obs, if (has_grid) grid_n else 0L)),
           file.path(dir, "meta.txt"))
cat(sprintf("wrote %d simulations (dim %d) and observation sets %s to %s\n",
            n, dim, paste(n_obs, collapse = ", "), dir))
if (has_grid) {
  cat(sprintf("wrote a %d^%d grid around the true parameter to %s\n",
              grid_n, dim, file.path(dir, "grid_theta.csv")))
}
