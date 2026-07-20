#' Posterior diagnostics
#'
#' Tools to check whether a trained posterior is trustworthy:
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
#' the true parameter within posterior samples conditioned on that data. If the
#' posterior is well calibrated, the ranks are uniformly distributed.
#'
#' @param fit An `nsbi_npe` fit (amortized posterior).
#' @param simulator The simulator used for inference.
#' @param prior The prior used for inference (defaults to `fit$prior`).
#' @param n_sbc Number of SBC trials (fresh (theta, x) pairs).
#' @param n_posterior_samples Posterior draws per trial (rank resolution).
#' @param seed Optional seed.
#' @return An object of class `nsbi_sbc` with the rank matrix and a per-parameter
#'   uniformity test.
#' @export
sbc <- function(fit, simulator, prior = fit$prior, n_sbc = 200L,
                n_posterior_samples = 1000L, seed = NULL) {
  stopifnot(inherits(fit, "nsbi_npe"))
  if (!is.null(seed)) set.seed(seed)
  d <- fit$dim_theta
  ranks <- matrix(NA_real_, nrow = n_sbc, ncol = d)
  theta_true <- sample_prior(prior, n_sbc)
  x_all <- as_theta_matrix(simulator(theta_true), fit$dim_x)
  for (i in seq_len(n_sbc)) {
    post <- posterior(fit, x_obs = x_all[i, ])
    draws <- sample.nsbi_posterior(post, n = n_posterior_samples)
    ranks[i, ] <- colSums(sweep(draws, 2, theta_true[i, ], `<`))
  }
  # per-parameter uniformity via chi-square on binned ranks
  L <- n_posterior_samples
  pvals <- apply(ranks, 2, function(r) {
    nb <- min(20L, L + 1L)
    br <- cut(r, breaks = seq(0, L, length.out = nb + 1L),
              include.lowest = TRUE)
    tab <- table(br)
    stats::chisq.test(tab)$p.value
  })
  structure(
    list(ranks = ranks, n_posterior_samples = L, n_sbc = n_sbc,
         uniformity_pvalue = pvals),
    class = "nsbi_sbc"
  )
}

#' @export
print.nsbi_sbc <- function(x, ...) {
  cat(sprintf("<nsbi_sbc> %d trials, %d posterior samples each\n",
              x$n_sbc, x$n_posterior_samples))
  cat("  per-parameter uniformity p-values (large = calibrated):\n")
  cat("   ", paste(sprintf("%.3f", x$uniformity_pvalue), collapse = "  "), "\n")
  invisible(x)
}

#' Expected coverage of central credible intervals
#'
#' Uses the SBC ranks to compare nominal credible levels with the empirical
#' fraction of trials in which the true parameter falls inside the corresponding
#' central interval. Well-calibrated posteriors lie on the diagonal.
#'
#' @param sbc_result An `nsbi_sbc` object from [sbc()].
#' @param levels Nominal credibility levels to evaluate.
#' @return A data frame with `nominal` and per-parameter empirical coverage.
#' @export
expected_coverage <- function(sbc_result, levels = seq(0.05, 0.95, by = 0.05)) {
  stopifnot(inherits(sbc_result, "nsbi_sbc"))
  L <- sbc_result$n_posterior_samples
  u <- sbc_result$ranks / L  # approx posterior CDF at truth ~ Uniform(0,1)
  emp <- sapply(levels, function(a) {
    lo <- (1 - a) / 2
    hi <- 1 - lo
    colMeans(u > lo & u < hi)
  })
  emp <- t(emp)
  colnames(emp) <- paste0("param", seq_len(ncol(emp)))
  data.frame(nominal = levels, emp, row.names = NULL, check.names = FALSE)
}

#' TARP expected coverage
#'
#' Tests of Accuracy with Random Points (Lemos et al. 2023). For each trial a
#' true parameter is drawn from the prior, data are simulated, and posterior
#' samples are drawn conditioned on those data. Given a random reference point,
#' the fraction of posterior samples closer to the reference than the truth is
#' the credibility level of the smallest distance-based credible region that
#' contains the truth. For a calibrated posterior these fractions are uniform,
#' so the expected coverage probability (ECP) at credibility level alpha equals
#' alpha.
#'
#' Unlike [sbc()], which ranks each parameter marginally, TARP is a *joint*
#' test: it can detect posteriors whose marginals are calibrated but whose
#' correlation structure is wrong. Distances are computed after z-scoring each
#' parameter (using the spread of the true draws), so parameters on different
#' scales contribute comparably.
#'
#' @param fit An `nsbi_npe` fit (amortized posterior).
#' @param simulator The simulator used for inference.
#' @param prior The prior used for inference (defaults to `fit$prior`).
#' @param n_tarp Number of TARP trials (fresh (theta, x) pairs).
#' @param n_posterior_samples Posterior draws per trial.
#' @param references How to draw reference points: `"uniform"` (default, uniform
#'   over the hyper-rectangle spanned by the true parameter draws, as in the
#'   paper) or `"prior"` (draws from the prior).
#' @param seed Optional seed.
#' @return An object of class `nsbi_tarp` with the per-trial coverage values
#'   and the ECP curve. Plot it with [plot_tarp()].
#' @references Lemos, Coogan, Hezaveh & Perreault-Levasseur (2023),
#'   "Sampling-based accuracy testing of posterior estimators for general
#'   inference", ICML. \doi{10.48550/arXiv.2302.03026}
#' @export
tarp <- function(fit, simulator, prior = fit$prior, n_tarp = 200L,
                 n_posterior_samples = 1000L,
                 references = c("uniform", "prior"), seed = NULL) {
  stopifnot(inherits(fit, "nsbi_npe"))
  references <- match.arg(references)
  if (!is.null(seed)) set.seed(seed)
  d <- fit$dim_theta

  theta_true <- sample_prior(prior, n_tarp)
  x_all <- as_theta_matrix(simulator(theta_true), fit$dim_x)

  # z-score all distances by the spread of the true draws so no single
  # parameter dominates
  std <- fit_standardizer(theta_true)
  theta_z <- apply_standardizer(std, theta_true)

  ref <- switch(references,
    uniform = {
      lo <- apply(theta_z, 2, min)
      hi <- apply(theta_z, 2, max)
      matrix(stats::runif(n_tarp * d, rep(lo, each = n_tarp),
                          rep(hi, each = n_tarp)), ncol = d)
    },
    prior = apply_standardizer(std, sample_prior(prior, n_tarp))
  )

  f <- numeric(n_tarp)
  for (i in seq_len(n_tarp)) {
    post <- posterior(fit, x_obs = x_all[i, ])
    draws <- sample.nsbi_posterior(post, n = n_posterior_samples)
    draws_z <- apply_standardizer(std, draws)
    d_samples <- sqrt(rowSums(sweep(draws_z, 2, ref[i, ], `-`)^2))
    d_truth <- sqrt(sum((theta_z[i, ] - ref[i, ])^2))
    f[i] <- mean(d_samples < d_truth)
  }

  levels <- seq(0, 1, by = 0.05)
  ecp <- sapply(levels, function(a) mean(f < a))
  structure(
    list(coverage_values = f, levels = levels, ecp = ecp,
         n_tarp = n_tarp, n_posterior_samples = n_posterior_samples,
         references = references),
    class = "nsbi_tarp"
  )
}

#' @export
print.nsbi_tarp <- function(x, ...) {
  cat(sprintf("<nsbi_tarp> %d trials, %d posterior samples each\n",
              x$n_tarp, x$n_posterior_samples))
  dev <- max(abs(x$ecp - x$levels))
  cat(sprintf("  max |ECP - nominal|: %.3f (0 = perfectly calibrated)\n", dev))
  cat("  plot with plot_tarp()\n")
  invisible(x)
}

#' Classifier two-sample test (C2ST)
#'
#' Trains a logistic-regression classifier to distinguish samples in `x` from
#' samples in `y` using cross-validation. A test accuracy near 0.5 means the two
#' sample sets are indistinguishable (good); near 1.0 means they differ. This is
#' the standard SBI metric for comparing an estimated posterior to a reference
#' (e.g. an analytic posterior or long-run MCMC draws).
#'
#' @param x,y Matrices of samples (rows = draws, cols = dimensions).
#' @param n_folds Number of cross-validation folds.
#' @param seed Optional seed.
#' @return A list with mean CV accuracy and per-fold accuracies.
#' @export
c2st <- function(x, y, n_folds = 5L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x <- as_theta_matrix(x)
  y <- as_theta_matrix(y)
  # standardize jointly for a fair, scale-free classifier
  data <- rbind(x, y)
  std <- fit_standardizer(data)
  data <- apply_standardizer(std, data)
  label <- c(rep(1, nrow(x)), rep(0, nrow(y)))
  n <- nrow(data)
  fold <- sample(rep_len(seq_len(n_folds), n))
  df <- data.frame(y = label, data)
  accs <- numeric(n_folds)
  for (k in seq_len(n_folds)) {
    tr <- fold != k
    te <- !tr
    fit <- suppressWarnings(
      stats::glm(y ~ ., data = df[tr, , drop = FALSE], family = stats::binomial())
    )
    pred <- stats::predict(fit, newdata = df[te, , drop = FALSE],
                           type = "response") > 0.5
    accs[k] <- mean(pred == (label[te] == 1))
  }
  list(accuracy = mean(accs), fold_accuracy = accs,
       interpretation = if (mean(accs) < 0.55)
         "indistinguishable (good)" else "distinguishable")
}

#' Posterior predictive draws
#'
#' Samples parameters from the posterior and pushes them back through the
#' simulator, giving predictive data to compare against the observation.
#'
#' @param post An `nsbi_posterior` object.
#' @param simulator The simulator.
#' @param n Number of predictive draws.
#' @param x Observation to condition on (defaults to `x_obs`).
#' @return An `n x d` matrix of simulated data from posterior parameter draws.
#' @export
posterior_predictive <- function(post, simulator, n = 1000L, x = NULL) {
  stopifnot(inherits(post, "nsbi_posterior"))
  theta <- sample.nsbi_posterior(post, n = n, obs = x)
  as_theta_matrix(simulator(theta))
}
