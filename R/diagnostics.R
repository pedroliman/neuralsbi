#' Posterior diagnostics
#'
#' Tools to check whether a trained posterior is trustworthy, mirroring the
#' validation utilities in the Python `sbi` package:
#'
#' * [sbc()] -- Simulation-Based Calibration rank statistics
#' * [expected_coverage()] -- nominal vs. empirical credible-interval coverage
#' * [c2st()] -- classifier two-sample test between two sample sets
#' * [posterior_predictive()] -- draw data from the fitted posterior
#'
#' @name diagnostics
NULL

#' Simulation-Based Calibration (SBC)
#'
#' Repeatedly draws a "true" parameter from the prior, simulates data, and ranks
#' the true parameter within posterior samples conditioned on that data. A
#' well-calibrated posterior gives uniformly distributed ranks.
#'
#' @param fit A trained [NPE] fit (amortized posterior).
#' @param simulator The simulator used for inference.
#' @param prior The prior used for inference (defaults to `fit@prior`).
#' @param n_sbc Number of SBC trials.
#' @param n_posterior_samples Posterior draws per trial.
#' @param seed Optional seed.
#' @return An [SBCResult] object.
#' @export
sbc <- function(fit, simulator, prior = fit@prior, n_sbc = 200L,
                n_posterior_samples = 1000L, seed = NULL) {
  if (!S7_inherits(fit, NPE)) stop("`fit` must be an NPE object.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  d <- fit@prior@n_dim
  ranks <- matrix(NA_real_, nrow = n_sbc, ncol = d)
  theta_true <- sample_prior(prior, n_sbc)
  x_all <- as_theta_matrix(simulator(theta_true), ncol(fit@x))
  for (i in seq_len(n_sbc)) {
    post <- posterior(fit, x_obs = x_all[i, ])
    draws <- sample(post, n = n_posterior_samples)
    ranks[i, ] <- colSums(sweep(draws, 2, theta_true[i, ], `<`))
  }
  L <- n_posterior_samples
  pvals <- apply(ranks, 2, function(r) {
    nb <- min(20L, L + 1L)
    tab <- table(cut(r, breaks = seq(0, L, length.out = nb + 1L),
                     include.lowest = TRUE))
    stats::chisq.test(tab)$p.value
  })
  SBCResult(ranks = ranks, n_posterior_samples = as.integer(L),
            n_sbc = as.integer(n_sbc), uniformity_pvalue = pvals)
}

method(print, SBCResult) <- function(x, ...) {
  cat(sprintf("<SBCResult> %d trials, %d posterior samples each\n",
              x@n_sbc, x@n_posterior_samples))
  cat("  per-parameter uniformity p-values (large = calibrated):\n")
  cat("   ", paste(sprintf("%.3f", x@uniformity_pvalue), collapse = "  "), "\n")
  invisible(x)
}

#' Expected coverage of central credible intervals
#'
#' Uses SBC ranks to compare nominal credible levels with the empirical fraction
#' of trials in which the true parameter lies inside the central interval.
#'
#' @param sbc_result An [SBCResult] from [sbc()].
#' @param levels Nominal credibility levels.
#' @return A data frame with `nominal` and per-parameter empirical coverage.
#' @export
expected_coverage <- function(sbc_result, levels = seq(0.05, 0.95, by = 0.05)) {
  if (!S7_inherits(sbc_result, SBCResult)) {
    stop("`sbc_result` must be an SBCResult.", call. = FALSE)
  }
  L <- sbc_result@n_posterior_samples
  u <- sbc_result@ranks / L
  emp <- t(sapply(levels, function(a) {
    lo <- (1 - a) / 2
    colMeans(u > lo & u < 1 - lo)
  }))
  colnames(emp) <- paste0("param", seq_len(ncol(emp)))
  data.frame(nominal = levels, emp, row.names = NULL, check.names = FALSE)
}

#' Classifier two-sample test (C2ST)
#'
#' Trains a cross-validated logistic-regression classifier to distinguish
#' samples in `x` from samples in `y`. Test accuracy near 0.5 means the two sets
#' are indistinguishable (good); near 1.0 means they differ. The standard SBI
#' metric for comparing an estimated posterior to a reference.
#'
#' @param x,y Matrices of samples (rows = draws).
#' @param n_folds Number of cross-validation folds.
#' @param seed Optional seed.
#' @return A list with mean CV accuracy, per-fold accuracies, interpretation.
#' @export
c2st <- function(x, y, n_folds = 5L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x <- as_theta_matrix(x); y <- as_theta_matrix(y)
  data <- rbind(x, y)
  std <- fit_standardizer(data)
  data <- apply_standardizer(std, data)
  label <- c(rep(1, nrow(x)), rep(0, nrow(y)))
  n <- nrow(data)
  fold <- base::sample(rep_len(seq_len(n_folds), n))
  df <- data.frame(y = label, data)
  accs <- numeric(n_folds)
  for (k in seq_len(n_folds)) {
    tr <- fold != k; te <- !tr
    fit <- suppressWarnings(
      stats::glm(y ~ ., data = df[tr, , drop = FALSE], family = stats::binomial()))
    pred <- stats::predict(fit, newdata = df[te, , drop = FALSE],
                           type = "response") > 0.5
    accs[k] <- mean(pred == (label[te] == 1))
  }
  list(accuracy = mean(accs), fold_accuracy = accs,
       interpretation = if (mean(accs) < 0.55) "indistinguishable (good)" else
         "distinguishable")
}

#' Posterior predictive draws
#'
#' Samples parameters from the posterior and pushes them through the simulator.
#'
#' @param post A [DirectPosterior] object.
#' @param simulator The simulator.
#' @param n Number of predictive draws.
#' @param x Observation to condition on (defaults to `x_obs`).
#' @return An `n x d` matrix of simulated data from posterior parameter draws.
#' @export
posterior_predictive <- function(post, simulator, n = 1000L, x = NULL) {
  theta <- sample(post, n = n, obs = x)
  as_theta_matrix(simulator(theta))
}
