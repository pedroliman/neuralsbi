#!/usr/bin/env Rscript
# Score both NLE implementations against the conjugate posterior, and each
# other.
#
# The reference is exact for every observation-set size: with prior N(0, v0 I)
# and likelihood N(theta, s2 I), n independent observations give a Gaussian
# posterior with variance 1 / (1/v0 + n/s2) and mean that variance times
# sum(x) / s2. So `c2st_vs_reference` is an accuracy number, not an agreement
# number, and the two implementations can be wrong in different directions
# without either being the yardstick.
#
# Usage: Rscript 08_compare_nle.R --estimator maf
suppressMessages(library(neuralsbi))

args <- as.list(commandArgs(trailingOnly = TRUE))
opt <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i)) default else args[[i + 1L]]
}
estimator <- opt("--estimator", "maf")
prior_var <- as.numeric(opt("--prior_var", "0.1"))
noise_var <- as.numeric(opt("--noise_var", "0.1"))

data_dir <- file.path("data", "gaussian_linear_iid")
res_dir <- file.path("results", "gaussian_linear_iid")
meta <- as.integer(readLines(file.path(data_dir, "meta.txt")))
n_obs <- meta[-1]

reference_iid <- function(x_obs, n) {
  post_var <- 1 / (1 / prior_var + nrow(x_obs) / noise_var)
  post_mean <- post_var * colSums(x_obs) / noise_var
  matrix(stats::rnorm(n * length(post_mean),
                      mean = rep(post_mean, each = n), sd = sqrt(post_var)),
         nrow = n)
}

rows <- list()
for (k in n_obs) {
  f_ours <- sprintf("%s/neuralsbi_%s_n%d.csv", res_dir, estimator, k)
  f_sbi <- sprintf("%s/sbi_%s_n%d.csv", res_dir, estimator, k)
  if (!file.exists(f_ours)) next
  x_obs <- as.matrix(read.csv(sprintf("%s/x_obs_n%d.csv", data_dir, k),
                              header = FALSE))
  ours <- as.matrix(read.csv(f_ours, header = FALSE))
  ref <- reference_iid(x_obs, nrow(ours))

  row <- data.frame(
    n_obs = k,
    c2st_ours_vs_ref = c2st(ours, ref, seed = k)$accuracy,
    ours_max_mean_err = max(abs(colMeans(ours) - colMeans(ref))),
    ours_sd_ratio = mean(apply(ours, 2, sd) / apply(ref, 2, sd))
  )
  if (file.exists(f_sbi)) {
    theirs <- as.matrix(read.csv(f_sbi, header = FALSE))
    ref_s <- reference_iid(x_obs, nrow(theirs))
    row$c2st_sbi_vs_ref <- c2st(theirs, ref_s, seed = k)$accuracy
    row$sbi_max_mean_err <- max(abs(colMeans(theirs) - colMeans(ref_s)))
    row$sbi_sd_ratio <- mean(apply(theirs, 2, sd) / apply(ref_s, 2, sd))
    row$c2st_ours_vs_sbi <- c2st(ours, theirs, seed = k)$accuracy
  }
  rows[[length(rows) + 1L]] <- row
}
tab <- do.call(rbind, rows)
print(tab, digits = 3)

out <- file.path(res_dir, sprintf("comparison_nle_%s.csv", estimator))
write.csv(tab, out, row.names = FALSE)
cat("wrote", out, "\n")
cat(sprintf("PASS criterion: c2st_ours_vs_ref <= 0.60 -> %s\n",
            if (all(tab$c2st_ours_vs_ref <= 0.60)) "PASS" else "FAIL"))
