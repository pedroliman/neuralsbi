#!/usr/bin/env Rscript
# Score both NRE implementations against the conjugate posterior, each other,
# and (when a grid was generated) the analytic log ratio.
#
# The reference is exact for every observation-set size, by the same formula
# 08_compare_nle.R uses: with prior N(0, v0 I) and likelihood N(theta, s2 I),
# n independent observations give a Gaussian posterior with variance
# 1 / (1/v0 + n/s2) and mean that variance times sum(x) / s2. c2st_*_vs_ref is
# therefore an accuracy number, not an agreement number.
#
# The grid piece is the half of #146 that was missing before: correlation and
# centred RMSE against the analytic log ratio say the two classifiers learned
# the same shape, but only c2st_ours_vs_sbi says the two posteriors are the
# same distribution, which is the actual M3 bar.
#
# Usage: Rscript 12_compare_nre.R --classifier resnet
suppressMessages(library(neuralsbi))

args <- as.list(commandArgs(trailingOnly = TRUE))
opt <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i)) default else args[[i + 1L]]
}
classifier <- opt("--classifier", "resnet")
prior_var <- as.numeric(opt("--prior_var", "0.1"))
noise_var <- as.numeric(opt("--noise_var", "0.1"))

data_dir <- file.path("data", "gaussian_linear_iid_nre")
res_dir <- file.path("results", "gaussian_linear_iid_nre")
meta <- as.integer(readLines(file.path(data_dir, "meta.txt")))
n_obs <- meta[2:(length(meta) - 1)]

reference_iid <- function(x_obs, n) {
  post_var <- 1 / (1 / prior_var + nrow(x_obs) / noise_var)
  post_mean <- post_var * colSums(x_obs) / noise_var
  matrix(stats::rnorm(n * length(post_mean),
                      mean = rep(post_mean, each = n), sd = sqrt(post_var)),
         nrow = n)
}

rows <- list()
for (k in n_obs) {
  f_ours <- sprintf("%s/neuralsbi_%s_n%d.csv", res_dir, classifier, k)
  f_sbi <- sprintf("%s/sbi_%s_n%d.csv", res_dir, classifier, k)
  if (!file.exists(f_ours)) next
  x_obs <- as.matrix(read.csv(sprintf("%s/x_obs_n%d.csv", data_dir, k),
                              header = FALSE))
  ours <- as.matrix(read.csv(f_ours, header = FALSE))
  ref <- reference_iid(x_obs, nrow(ours))

  row <- data.frame(
    n_obs = k,
    c2st_ours_vs_ref = c2st(ref, ours, seed = k)$accuracy,
    ours_max_mean_err = max(abs(colMeans(ours) - colMeans(ref))),
    ours_sd_ratio = mean(apply(ours, 2, sd) / apply(ref, 2, sd))
  )
  if (file.exists(f_sbi)) {
    theirs <- as.matrix(read.csv(f_sbi, header = FALSE))
    ref_s <- reference_iid(x_obs, nrow(theirs))
    row$c2st_sbi_vs_ref <- c2st(ref_s, theirs, seed = k)$accuracy
    row$sbi_max_mean_err <- max(abs(colMeans(theirs) - colMeans(ref_s)))
    row$sbi_sd_ratio <- mean(apply(theirs, 2, sd) / apply(ref_s, 2, sd))
    row$c2st_ours_vs_sbi <- c2st(ours, theirs, seed = k)$accuracy
  }
  rows[[length(rows) + 1L]] <- row
}
tab <- do.call(rbind, rows)
print(tab, digits = 3)

out <- file.path(res_dir, sprintf("comparison_nre_%s.csv", classifier))
write.csv(tab, out, row.names = FALSE)
cat("wrote", out, "\n")
cat(sprintf("PASS criterion: c2st_ours_vs_ref <= 0.60 -> %s\n",
            if (all(tab$c2st_ours_vs_ref <= 0.60)) "PASS" else "FAIL"))
if ("c2st_ours_vs_sbi" %in% names(tab)) {
  cat(sprintf("PASS criterion: c2st_ours_vs_sbi <= 0.60 -> %s\n",
              if (all(tab$c2st_ours_vs_sbi <= 0.60)) "PASS" else "FAIL"))
}

# Secondary metric (#146): correlation and centred RMSE of the learned ratio
# against the analytic one, at the first observation. Both implementations'
# ratios are identified only up to an additive constant (see ?nre), so each
# series is centred on its own mean before comparing. analytic_log_ratio drops
# every term that does not depend on theta, which is all a ratio identifies
# anyway -- see task_gaussian_linear()'s simulator for the likelihood this
# comes from. The grid file exists only when 09_generate_data_nre.R judged the
# dimension low enough to grid at all.
grid_path <- file.path(data_dir, "grid_theta.csv")
if (file.exists(grid_path)) {
  grid_theta <- as.matrix(read.csv(grid_path, header = FALSE))
  x_obs_1 <- as.matrix(read.csv(file.path(data_dir, "x_obs_n1.csv"),
                                header = FALSE))
  analytic <- -0.5 * rowSums(sweep(grid_theta, 2, as.numeric(x_obs_1), "-")^2) /
    noise_var
  analytic_c <- analytic - mean(analytic)

  grid_rows <- list()
  for (side in c("neuralsbi", "sbi")) {
    f <- file.path(res_dir, sprintf("grid_log_ratio_%s_%s.csv", side, classifier))
    if (!file.exists(f)) next
    est <- as.numeric(read.csv(f, header = FALSE)[[1]])
    est_c <- est - mean(est)
    grid_rows[[side]] <- data.frame(
      source = side,
      corr = stats::cor(est_c, analytic_c),
      centred_rmse = sqrt(mean((est_c - analytic_c)^2))
    )
  }
  if (length(grid_rows)) {
    grid_tab <- do.call(rbind, grid_rows)
    cat("\nanalytic-log-ratio grid comparison:\n")
    print(grid_tab, digits = 3)
    grid_out <- file.path(res_dir,
                          sprintf("comparison_nre_grid_%s.csv", classifier))
    write.csv(grid_tab, grid_out, row.names = FALSE)
    cat("wrote", grid_out, "\n")
  }
}
