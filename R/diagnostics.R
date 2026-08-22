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

#' Accept any kind of fit, and say which ones exist when none matches
#' @keywords internal
check_inference_fit <- function(fit) {
  if (inherits(fit, "nsbi_npe") || inherits(fit, "nsbi_nle") ||
      inherits(fit, "nsbi_nre")) {
    return(invisible(TRUE))
  }
  stop("Expected a fit from npe(), nle() or nre(), not an object of class ",
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

#' Shared simulate-and-drop preamble for sbc() and tarp()
#'
#' Both diagnostics draw "true" parameters from `prior`, simulate data for
#' them, and drop trials whose simulation failed. The prior-width check lives
#' here rather than in each caller, so `sbc()` and `tarp()` get it in one
#' place: `prior` exists to be overridden (see their docs), and a prior of the
#' wrong width would otherwise reach `sweep()`/z-scoring downstream and
#' silently recycle against the fit's width.
#' @param fit An `nsbi_npe`, `nsbi_nle` or `nsbi_nre` fit.
#' @param simulator The simulator used for inference.
#' @param prior The prior to draw true parameters from.
#' @param n Number of trials to draw.
#' @param sim_args Named list of extra simulator arguments.
#' @param what Label for [drop_failed_sims()]'s warning (e.g. `"SBC trials"`).
#' @return A list with `theta` and `x` (the surviving trials), `n` (their
#'   count, i.e. `nrow(theta)` after dropping), and `n_dropped`.
#' @keywords internal
sbc_draws <- function(fit, simulator, prior, n, sim_args, what) {
  check_inference_fit(fit)
  check_prior(prior, dim = fit$dim_theta)
  theta_true <- sample_prior(prior, n)
  x_all <- run_simulator(simulator, theta_true, sim_args = sim_args,
                         d = fit$dim_x)
  kept <- drop_failed_sims(theta_true, x_all, what = what)
  list(theta = kept$theta, x = kept$x, n = nrow(kept$theta),
       n_dropped = kept$n_dropped)
}

#' Shared per-trial posterior-draw loop for sbc() and tarp()
#'
#' Draws a posterior for each row of `x_all`, insists on the full
#' `n_posterior_samples` count via [diagnostic_draws()], and hands the result
#' to `f(draws, i)` for whichever per-trial metric the caller wants: `sbc()`
#' returns a rank row, `tarp()` a scalar coverage value.
#'
#' The denominator a caller bins or normalizes against (`sbc()`'s rank scale,
#' via the p-value bins and [expected_coverage()]) is read back from how many
#' draws a trial actually returned, not trusted from the `n_posterior_samples`
#' argument -- so if a future change ever lets a short draw through instead of
#' erroring, the diagnostic is still scored on the scale it was actually drawn
#' on rather than the one it was asked for.
#'
#' `finally`, not `on.exit()`: this block runs inside
#' `progressr::with_progress()`, so `on.exit()` would attach to the promise's
#' forcing frame and close the bar before the loop starts (see
#' `run_simulator()` in R/parallel.R).
#' @param fit An `nsbi_npe`, `nsbi_nle` or `nsbi_nre` fit.
#' @param x_all Matrix of simulated data, one row (one trial) per observation.
#' @param n_posterior_samples Posterior draws requested per trial.
#' @param label Progress-bar label.
#' @param f `function(draws, i)`, the per-trial metric.
#' @param ... Passed to [posterior()].
#' @return A list with `results` (a list of length `nrow(x_all)`, one `f()`
#'   return value per trial) and `n_posterior_samples` (the draw count trials
#'   were actually scored against).
#' @keywords internal
for_each_trial <- function(fit, x_all, n_posterior_samples, label, f, ...) {
  n <- nrow(x_all)
  results <- vector("list", n)
  n_drawn <- n_posterior_samples
  with_nsbi_progress({
    p <- nsbi_progressor(steps = n, label = label)
    tryCatch({
      for (i in seq_len(n)) {
        post <- posterior(fit, x_obs = x_all[i, ], ...)
        draws <- diagnostic_draws(post, n_posterior_samples, i)
        n_drawn <- nrow(draws)
        results[[i]] <- f(draws, i)
        p(1)
      }
    }, finally = p(0, done = TRUE))
  })
  list(results = results, n_posterior_samples = n_drawn)
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
#' @param fit An `nsbi_npe` fit from [npe()], an `nsbi_nle` fit from [nle()],
#'   or an `nsbi_nre` fit from [nre()]. With an NLE or NRE fit every trial is a
#'   separate MCMC run, so start with a small `n_sbc` and raise it once the
#'   cost is known.
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
#'   (`n_chains`, `warmup`, `thin`, `sampler`) reach an NLE or NRE fit.
#' @details The `n_sbc` simulations run across `future` workers when a plan is
#'   set (see [nsbi_parallel]); the ranking loop that follows calls the trained
#'   network and always runs locally.
#' @return An object of class `nsbi_sbc` with the rank matrix and a per-parameter
#'   uniformity test.
#' @export
sbc <- function(fit, simulator, prior = fit$prior, n_sbc = 200L,
                n_posterior_samples = 1000L, sim_args = list(),
                seed = NULL, ...) {
  n_sbc <- check_count(n_sbc, "n_sbc")
  n_posterior_samples <- check_count(n_posterior_samples,
                                     "n_posterior_samples")
  if (!is.null(seed)) set.seed(seed)
  prep <- sbc_draws(fit, simulator, prior, n_sbc, sim_args, what = "SBC trials")
  theta_true <- prep$theta
  n_sbc <- prep$n

  trial <- for_each_trial(fit, prep$x, n_posterior_samples, "Ranking",
    function(draws, i) colSums(sweep(draws, 2, theta_true[i, ], `<`)), ...)
  ranks <- do.call(rbind, trial$results) %||%
    matrix(NA_real_, nrow = 0L, ncol = fit$dim_theta)
  if (is.null(colnames(ranks)) && !is.null(fit$param_names)) {
    colnames(ranks) <- fit$param_names
  }
  # per-parameter uniformity via chi-square on binned ranks. Binned against
  # what trials were actually scored against (see for_each_trial()), not the
  # requested n_posterior_samples -- these agree whenever every trial hits its
  # full draw count, which diagnostic_draws() currently enforces, but the rank
  # scale should follow the draws rather than the request either way.
  L <- trial$n_posterior_samples
  pvals <- apply(ranks, 2, function(r) {
    nb <- min(20L, L + 1L)
    br <- cut(r, breaks = seq(0, L, length.out = nb + 1L),
              include.lowest = TRUE)
    tab <- table(br)
    # Monte Carlo p-value: bypasses the asymptotic chi-squared approximation
    # entirely (and its low-expected-count warning), rather than working
    # around a case where the approximation happens to hold. Also valid at
    # the small n_sbc test fixtures use, where the asymptotic version warns
    # routinely.
    stats::chisq.test(tab, simulate.p.value = TRUE)$p.value
  })
  names(pvals) <- colnames(ranks)
  structure(
    list(ranks = ranks, n_posterior_samples = L, n_sbc = n_sbc,
         n_dropped = prep$n_dropped, uniformity_pvalue = pvals),
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
  emp <- vapply(levels, function(a) {
    lo <- (1 - a) / 2
    hi <- 1 - lo
    colMeans(u > lo & u < hi)
  }, numeric(ncol(u)))
  # vapply drops the params x levels matrix to a plain length(levels) vector
  # when ncol(u) == 1 (a single-parameter fit), since each call's own return
  # value is then also length 1. Put the dropped dimension back before
  # transposing, or the lone parameter's per-level coverage gets read back out
  # column-major as one row (i.e. one fake parameter) per nominal level below,
  # each holding the same three numbers repeated across every level.
  if (is.null(dim(emp))) dim(emp) <- c(ncol(u), length(levels))
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
#' @param fit An `nsbi_npe` fit from [npe()], an `nsbi_nle` fit from [nle()],
#'   or an `nsbi_nre` fit from [nre()]. With an NLE or NRE fit every trial is a
#'   separate MCMC run, so start with a small `n_tarp` and raise it once the
#'   cost is known.
#' @param simulator The simulator used for inference; called once per trial
#'   (see [nsbi_simulator]).
#' @param prior The prior to draw the true parameters from, and the reference
#'   points when `references = "prior"` (defaults to `fit$prior`). As in
#'   [sbc()], the coverage claim is about the prior the fit was trained on, so
#'   the default is what tests this fit; another prior tests how the fit does on
#'   parameters it was not calibrated against. It must cover the same
#'   parameters as the fit.
#' @param ... Passed to [posterior()], which is how the MCMC controls
#'   (`n_chains`, `warmup`, `thin`, `sampler`) reach an NLE or NRE fit.
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
  references <- match.arg(references)
  n_tarp <- check_count(n_tarp, "n_tarp")
  n_posterior_samples <- check_count(n_posterior_samples,
                                     "n_posterior_samples")
  if (!is.null(seed)) set.seed(seed)
  prep <- sbc_draws(fit, simulator, prior, n_tarp, sim_args, what = "TARP trials")
  theta_true <- prep$theta
  n_tarp <- prep$n
  d <- fit$dim_theta

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

  trial <- for_each_trial(fit, prep$x, n_posterior_samples, "Coverage",
    function(draws, i) {
      draws_z <- apply_standardizer(std, draws)
      d_samples <- sqrt(rowSums(sweep(draws_z, 2, ref[i, ], `-`)^2))
      d_truth <- sqrt(sum((theta_z[i, ] - ref[i, ])^2))
      mean(d_samples < d_truth)
    }, ...)
  f <- unlist(trial$results) %||% numeric(0)

  levels <- seq(0, 1, by = 0.05)
  ecp <- sapply(levels, function(a) mean(f < a))
  structure(
    list(coverage_values = f, levels = levels, ecp = ecp,
         n_tarp = n_tarp, n_dropped = prep$n_dropped,
         n_posterior_samples = trial$n_posterior_samples,
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
  if (nrow(theta) == 0L) {
    rate <- attr(theta, "acceptance_rate") %||% 0
    stop(sprintf(paste0(
      "posterior_predictive(): the posterior returned no draws inside the ",
      "prior support at this observation (acceptance rate %.2f), so there ",
      "is nothing to push through the simulator. This usually means the ",
      "observation is outside the range the fit was trained on."), rate),
      call. = FALSE)
  }
  pred <- run_simulator(simulator, theta, sim_args = sim_args,
                        label = "Predicting")
  pred <- drop_failed_sims(NULL, pred, what = "predictive draws")$x
  if (is.null(colnames(pred)) && !is.null(post$fit$x_names)) {
    colnames(pred) <- post$fit$x_names
  }
  pred
}
