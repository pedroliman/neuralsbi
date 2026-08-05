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

#' Accept either kind of fit, and say which ones exist when neither matches
#' @keywords internal
check_inference_fit <- function(fit) {
  if (inherits(fit, "nsbi_npe") || inherits(fit, "nsbi_nle")) {
    return(invisible(TRUE))
  }
  stop("Expected a fit from npe() or nle(), not an object of class ",
       paste(class(fit), collapse = "/"), ".", call. = FALSE)
}

#' Draw posterior samples for one diagnostic trial, insisting on the full count
#'
#' `sample.nsbi_posterior()` returns fewer rows than asked for when a bounded
#' prior and a leaky estimator leave rejection sampling short, and only warns.
#' The diagnostics cannot absorb that quietly. `sbc()` bins its ranks against
#' `n_posterior_samples`, so a trial that came back short is scored on a scale it
#' was never drawn on, and the ranks compress toward zero; rescaling that one
#' trial on its own would instead make it incomparable to the others. Either way
#' the run reports a miscalibrated posterior when the real cause is lost draws,
#' so stop and say so.
#' @keywords internal
diagnostic_draws <- function(post, n, trial) {
  draws <- sample(post, n = n)
  if (nrow(draws) < n) {
    stop(sprintf(
      paste0("Trial %d returned %d of %d posterior draws. The estimator is ",
             "leaking mass outside the prior support, and a short draw would ",
             "be scored against %d, biasing the diagnostic toward ",
             "'miscalibrated'. Train on more simulations, or lower ",
             "`n_posterior_samples`."),
      trial, nrow(draws), n, n), call. = FALSE)
  }
  draws
}

#' Simulation-Based Calibration (SBC)
#'
#' Repeatedly draws a "true" parameter from the prior, simulates data, and ranks
#' the true parameter within posterior samples conditioned on that data. If the
#' posterior is well calibrated, the ranks are uniformly distributed.
#'
#' A trial whose simulation returns non-finite output is dropped, which lowers
#' the effective `n_sbc`. A trial whose posterior comes back with fewer than
#' `n_posterior_samples` draws, which happens when a bounded prior and a leaky
#' estimator defeat rejection sampling, is an error: ranks are binned against
#' `n_posterior_samples`, so a short draw would be scored on a scale it was
#' never drawn on and would read as miscalibration.
#'
#' @param fit An `nsbi_npe` fit from [npe()], or an `nsbi_nle` fit from
#'   [nle()]. With an NLE fit every trial is a separate MCMC run, so start
#'   with a small `n_sbc` and raise it once the cost is known.
#' @param simulator The simulator used for inference; called once per trial
#'   (see [nsbi_simulator]).
#' @param prior The prior to draw the true parameters from (defaults to
#'   `fit$prior`). SBC is a test of the posterior against the prior it was
#'   trained on, so the default is the only choice that answers "is this fit
#'   calibrated". Overriding it changes the question to how the fit behaves on
#'   parameters drawn from somewhere else, which is a reasonable local check
#'   but is no longer SBC. It must cover the same parameters as the fit.
#' @param n_sbc Number of SBC trials (fresh (theta, x) pairs).
#' @param n_posterior_samples Posterior draws per trial (rank resolution).
#' @param sim_args Named list of extra arguments passed to every simulator
#'   call; see [nsbi_simulator].
#' @param seed Optional seed.
#' @param ... Passed to [posterior()], which is how the MCMC controls
#'   (`n_chains`, `warmup`, `thin`, `sampler`) reach an NLE fit.
#' @details The `n_sbc` simulations run across `future` workers when a plan is
#'   set (see [nsbi_parallel]); the ranking loop that follows calls the trained
#'   network and always runs locally.
#' @return An object of class `nsbi_sbc` with the rank matrix and a per-parameter
#'   uniformity test.
#' @export
sbc <- function(fit, simulator, prior = fit$prior, n_sbc = 200L,
                n_posterior_samples = 1000L, sim_args = list(),
                seed = NULL, ...) {
  check_inference_fit(fit)
  # The ranks are sized from the fit but the truths come from `prior`, so a
  # prior of the wrong width makes sweep() recycle one against the other and
  # the ranks come out of a comparison nobody asked for. Check it here, above
  # the simulation loop: n_sbc simulator calls is a long way to travel before
  # failing on an argument.
  check_prior(prior, dim = fit$dim_theta)
  n_sbc <- check_count(n_sbc, "n_sbc")
  n_posterior_samples <- check_count(n_posterior_samples,
                                     "n_posterior_samples")
  if (!is.null(seed)) set.seed(seed)
  d <- fit$dim_theta
  theta_true <- sample_prior(prior, n_sbc)
  x_all <- run_simulator(simulator, theta_true, sim_args = sim_args,
                         d = fit$dim_x)
  kept <- drop_failed_sims(theta_true, x_all, what = "SBC trials")
  theta_true <- kept$theta
  x_all <- kept$x
  n_sbc <- nrow(theta_true)
  ranks <- matrix(NA_real_, nrow = n_sbc, ncol = d)
  with_nsbi_progress({
    p <- nsbi_progressor(steps = n_sbc, label = "Ranking")
    tryCatch({
      for (i in seq_len(n_sbc)) {
        post <- posterior(fit, x_obs = x_all[i, ], ...)
        draws <- diagnostic_draws(post, n_posterior_samples, i)
        ranks[i, ] <- colSums(sweep(draws, 2, theta_true[i, ], `<`))
        p(1)
      }
    }, finally = p(0, done = TRUE))
  })
  if (is.null(colnames(ranks)) && !is.null(fit$param_names)) {
    colnames(ranks) <- fit$param_names
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
  names(pvals) <- colnames(ranks)
  structure(
    list(ranks = ranks, n_posterior_samples = L, n_sbc = n_sbc,
         n_dropped = kept$n_dropped, uniformity_pvalue = pvals),
    class = "nsbi_sbc"
  )
}

#' @export
print.nsbi_sbc <- function(x, ...) {
  cat(sprintf("<nsbi_sbc> %d trials, %d posterior samples each\n",
              x$n_sbc, x$n_posterior_samples))
  cat_dropped(x$n_dropped,
             what = "further trials dropped for non-finite simulator output")
  cat("  per-parameter uniformity p-values (large = calibrated):\n")
  if (!is.null(names(x$uniformity_pvalue))) {
    cat("   ", paste(names(x$uniformity_pvalue), sprintf("%.3f", x$uniformity_pvalue),
                     sep = "=", collapse = "  "), "\n")
  } else {
    cat("   ", paste(sprintf("%.3f", x$uniformity_pvalue), collapse = "  "), "\n")
  }
  invisible(x)
}

#' Expected coverage of central credible intervals
#'
#' Uses the SBC ranks to compare nominal credible levels with the empirical
#' fraction of trials in which the true parameter falls inside the corresponding
#' central interval. Well-calibrated posteriors lie on the diagonal.
#'
#' @param sbc_result An `nsbi_sbc` object from [sbc()].
#' @param levels Nominal credibility levels to evaluate, each strictly between
#'   0 and 1.
#' @return A data frame with `nominal` and per-parameter empirical coverage.
#' @export
expected_coverage <- function(sbc_result, levels = seq(0.05, 0.95, by = 0.05)) {
  stopifnot(inherits(sbc_result, "nsbi_sbc"))
  # A level outside (0, 1) is not a credible interval. It used to be scored
  # anyway: the interval comes out empty or covers the line, so the row reads
  # as coverage 0 or 1 and looks like a badly calibrated posterior.
  levels <- check_probs(levels, "levels")
  L <- sbc_result$n_posterior_samples
  u <- sbc_result$ranks / L  # approx posterior CDF at truth ~ Uniform(0,1)
  emp <- sapply(levels, function(a) {
    lo <- (1 - a) / 2
    hi <- 1 - lo
    colMeans(u > lo & u < hi)
  })
  emp <- t(emp)
  colnames(emp) <- colnames(sbc_result$ranks) %||% paste0("param", seq_len(ncol(emp)))
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
#' A trial whose simulation returns non-finite output is dropped, which lowers
#' the effective `n_tarp`. A trial whose posterior comes back with fewer than
#' `n_posterior_samples` draws is an error, for the same reason as in [sbc()]:
#' a trial scored on a different number of draws is not comparable to the rest.
#'
#' @param fit An `nsbi_npe` fit from [npe()], or an `nsbi_nle` fit from
#'   [nle()]. With an NLE fit every trial is a separate MCMC run, so start
#'   with a small `n_tarp` and raise it once the cost is known.
#' @param simulator The simulator used for inference; called once per trial
#'   (see [nsbi_simulator]).
#' @param prior The prior to draw the true parameters from, and the reference
#'   points when `references = "prior"` (defaults to `fit$prior`). As in
#'   [sbc()], the coverage claim is about the prior the fit was trained on, so
#'   the default is what tests this fit; another prior tests how the fit does on
#'   parameters it was not calibrated against. It must cover the same
#'   parameters as the fit.
#' @param ... Passed to [posterior()], which is how the MCMC controls
#'   (`n_chains`, `warmup`, `thin`, `sampler`) reach an NLE fit.
#' @param n_tarp Number of TARP trials (fresh (theta, x) pairs).
#' @param n_posterior_samples Posterior draws per trial.
#' @param references How to draw reference points: `"uniform"` (default, uniform
#'   over the hyper-rectangle spanned by the true parameter draws, as in the
#'   paper) or `"prior"` (draws from the prior).
#' @param sim_args Named list of extra arguments passed to every simulator
#'   call; see [nsbi_simulator].
#' @param seed Optional seed.
#' @return An object of class `nsbi_tarp` with the per-trial coverage values
#'   and the ECP curve. Plot it with [plot_tarp()].
#' @references Lemos, Coogan, Hezaveh & Perreault-Levasseur (2023),
#'   "Sampling-based accuracy testing of posterior estimators for general
#'   inference", ICML. \doi{10.48550/arXiv.2302.03026}
#' @export
tarp <- function(fit, simulator, prior = fit$prior, n_tarp = 200L,
                 n_posterior_samples = 1000L,
                 references = c("uniform", "prior"), sim_args = list(),
                 seed = NULL, ...) {
  check_inference_fit(fit)
  # Same reason as sbc(): the truths come from `prior` and everything they are
  # measured against is sized from the fit, here the z-scoring and the
  # distances. Above the simulation loop, so a wrong prior costs nothing.
  check_prior(prior, dim = fit$dim_theta)
  references <- match.arg(references)
  n_tarp <- check_count(n_tarp, "n_tarp")
  n_posterior_samples <- check_count(n_posterior_samples,
                                     "n_posterior_samples")
  if (!is.null(seed)) set.seed(seed)
  d <- fit$dim_theta

  theta_true <- sample_prior(prior, n_tarp)
  x_all <- run_simulator(simulator, theta_true, sim_args = sim_args,
                         d = fit$dim_x)
  kept <- drop_failed_sims(theta_true, x_all, what = "TARP trials")
  theta_true <- kept$theta
  x_all <- kept$x
  n_tarp <- nrow(theta_true)

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
  with_nsbi_progress({
    p <- nsbi_progressor(steps = n_tarp, label = "Coverage")
    tryCatch({
      for (i in seq_len(n_tarp)) {
        post <- posterior(fit, x_obs = x_all[i, ], ...)
        draws <- diagnostic_draws(post, n_posterior_samples, i)
        draws_z <- apply_standardizer(std, draws)
        d_samples <- sqrt(rowSums(sweep(draws_z, 2, ref[i, ], `-`)^2))
        d_truth <- sqrt(sum((theta_z[i, ] - ref[i, ])^2))
        f[i] <- mean(d_samples < d_truth)
        p(1)
      }
    }, finally = p(0, done = TRUE))
  })

  levels <- seq(0, 1, by = 0.05)
  ecp <- sapply(levels, function(a) mean(f < a))
  structure(
    list(coverage_values = f, levels = levels, ecp = ecp,
         n_tarp = n_tarp, n_dropped = kept$n_dropped,
         n_posterior_samples = n_posterior_samples,
         references = references),
    class = "nsbi_tarp"
  )
}

#' @export
print.nsbi_tarp <- function(x, ...) {
  cat(sprintf("<nsbi_tarp> %d trials, %d posterior samples each\n",
              x$n_tarp, x$n_posterior_samples))
  cat_dropped(x$n_dropped,
             what = "further trials dropped for non-finite simulator output")
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
#' The classifier is linear, which decides what the number can and cannot see.
#' It picks up a shift in location readily; it is close to blind to two sample
#' sets that share a mean and differ in spread or in the shape of their
#' dependence, because no hyperplane separates those. Python `sbi` uses an MLP,
#' so its C2ST sees more, and a 0.5 from here is the weaker claim of the two.
#' Read it alongside the moments rather than on its own.
#'
#' Unequal sample sizes are balanced by subsampling the larger set, because
#' accuracy against unbalanced classes is not a two-sample test: 8000 draws
#' against 2000 identical ones scores 0.8 for a classifier that has learned
#' nothing except to always answer with the bigger class.
#'
#' @param x,y Matrices of samples (rows = draws, cols = dimensions). The two
#'   must be the same width, since they are draws of the same quantity. Row
#'   counts need not match; the larger set is subsampled down to the smaller.
#' @param n_folds Number of cross-validation folds. At least 2, and fewer than
#'   the number of draws in the smaller sample set.
#' @param seed Optional seed.
#' @return A list with mean CV accuracy and per-fold accuracies.
#' @export
c2st <- function(x, y, n_folds = 5L, seed = NULL) {
  n_folds <- check_count(n_folds, "n_folds", min = 2L,
                         why = "since each fold is scored by a fit on the rest")
  if (!is.null(seed)) set.seed(seed)
  # A row is one draw here, so a bare vector is a column of 1-D draws rather
  # than check_matrix()'s single row. That is the pre-computed (theta, x) rule
  # in npe(), and it is why the type check and the reshape are separate calls.
  x <- as_theta_matrix(check_numeric(x, "x"))
  y <- as_theta_matrix(check_numeric(y, "y"))
  if (ncol(x) != ncol(y)) {
    # rbind() is where this lands otherwise, and it complains about the number
    # of columns of a matrix it was handed rather than about x or y.
    stop(sprintf(paste0("`x` has %s and `y` has %s. c2st() compares two sets ",
                        "of draws of the same quantity, so both must be the ",
                        "same width (rows = draws, columns = dimensions)."),
                 n_things(ncol(x), "column"), n_things(ncol(y), "column")),
         call. = FALSE)
  }
  n_each <- min(nrow(x), nrow(y))
  if (n_folds >= n_each) {
    # Folds are cut over both sample sets together, so n_folds up to n_each
    # still leaves each training fit some draws from each class. Past that the
    # test folds thin out and empty ones score NaN, which propagates to the
    # accuracy with nothing said about why.
    stop(sprintf(paste0("`n_folds` is %d, but the smaller sample set has only ",
                        "%s. It must be fewer, so that every fold has draws to ",
                        "test on."),
                 n_folds, n_things(n_each, "draw")),
         call. = FALSE)
  }
  # base::sample.int, not the package's sample() generic, which dispatches on
  # its first argument.
  if (nrow(x) > n_each) x <- x[base::sample.int(nrow(x), n_each), , drop = FALSE]
  if (nrow(y) > n_each) y <- y[base::sample.int(nrow(y), n_each), , drop = FALSE]
  # standardize jointly for a fair, scale-free classifier
  data <- rbind(x, y)
  std <- fit_standardizer(data)
  data <- apply_standardizer(std, data)
  label <- c(rep(1, nrow(x)), rep(0, nrow(y)))
  n <- nrow(data)
  # base::sample again: sample() is an S3 generic in this package.
  fold <- base::sample(rep_len(seq_len(n_folds), n))
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
#' simulator, giving predictive data to compare against the observation. A draw
#' whose simulation returns non-finite output is dropped, which lowers the
#' number of predictive draws returned.
#'
#' @param post An `nsbi_posterior` object.
#' @param simulator The simulator; called once per posterior draw (see
#'   [nsbi_simulator]).
#' @param n Number of predictive draws.
#' @param x Observation to condition on (defaults to `x_obs`).
#' @param sim_args Named list of extra arguments passed to every simulator
#'   call; see [nsbi_simulator].
#' @return An `n x d` matrix of simulated data from posterior parameter draws.
#' @export
posterior_predictive <- function(post, simulator, n = 1000L, x = NULL,
                                 sim_args = list()) {
  stopifnot(inherits(post, "nsbi_posterior"))
  theta <- sample(post, n = n, obs = x)
  pred <- run_simulator(simulator, theta, sim_args = sim_args,
                        label = "Predicting")
  pred <- drop_failed_sims(NULL, pred, what = "predictive draws")$x
  if (is.null(colnames(pred)) && !is.null(post$fit$x_names)) {
    colnames(pred) <- post$fit$x_names
  }
  pred
}
