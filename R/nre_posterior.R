#' Posterior from a neural ratio
#'
#' Bayes' rule turns an [nre()] fit into a posterior,
#' \eqn{p(\theta \mid x) \propto r_\phi(\theta, x)\,p(\theta)}, and as with
#' [nle()] the result has no closed form and no direct sampler. `posterior()`
#' on an `nsbi_nre` therefore returns an object that samples with MCMC, and
#' [sample()] on it runs a chain rather than a forward pass.
#'
#' There is one sampler here, the vectorized slice sampler of [nsbi_mcmc], and
#' no `sampler` argument to choose it with. [nle()]'s `"stan"` option works by
#' transpiling the fitted density into a Stan `functions` block; the residual
#' classifier behind a ratio estimator has no such export, so NUTS is not on
#' offer.
#'
#' @param fit An `nsbi_nre` object from [nre()].
#' @param x_obs Observation to condition on. Rows are treated as independent
#'   observations from the same parameter, and the log ratio sums over them.
#' @param n_chains Number of chains. The slice sampler evaluates every chain in
#'   one batched call per step, so more chains cost almost nothing.
#' @param warmup Steps discarded at the start of each chain.
#' @param thin Keep one draw in `thin`. See [posterior.nsbi_nle()] for why the
#'   default is 2.
#' @param init_strategy `"resample"` (default, and `sbi`'s) weights a pool of
#'   prior draws by the posterior density and resamples the starting points
#'   from it; `"proposal"` starts from plain prior draws.
#' @param seed Optional integer seed.
#' @param ... Further arguments to the sampler: `width`, `max_steps` and
#'   `n_pool`.
#'
#' @return An object of class `c("nsbi_nre_posterior", "nsbi_posterior")`.
#'
#' @examples
#' prior <- prior_uniform(c(mu = -3), c(mu = 3))
#' fit <- nre(prior, function(mu) c(y = rnorm(1, mu, 0.5)),
#'            n_simulations = 1000, classifier = "logistic")
#'
#' x_obs <- matrix(rnorm(50, mean = 1, sd = 0.5), ncol = 1)
#' post <- posterior(fit, x_obs, n_chains = 4, warmup = 50, thin = 2)
#' draws <- sample(post, 400)
#' colMeans(draws)
#' @export
posterior.nsbi_nre <- function(fit, x_obs = NULL, n_chains = 20L,
                               warmup = 200L, thin = 2L,
                               init_strategy = c("resample", "proposal"),
                               seed = NULL, ...) {
  check_fit_alive(fit)
  init_strategy <- match.arg(init_strategy)
  if (!is.null(x_obs)) {
    x_obs <- check_numeric(x_obs, "x_obs")
    check_finite(x_obs, "x_obs")
    x_obs <- as_theta_matrix(x_obs, fit$dim_x)
  }
  n_chains <- check_mcmc_count(n_chains, "n_chains", 2L,
                               "so convergence can be diagnosed")
  warmup <- check_mcmc_count(warmup, "warmup", 0L)
  thin <- check_mcmc_count(thin, "thin", 1L,
                           "since one draw in `thin` is kept")

  mcmc_posterior(fit, x_obs, "slice", n_chains, warmup, thin, init_strategy,
                 seed, list(...), "nsbi_nre_posterior")
}

#' Sample an NRE posterior with MCMC
#'
#' Identical in behavior to [sample.nsbi_nle_posterior()], including the draw
#' cache: asking for the same or fewer draws from the same observation returns
#' the cached run instead of starting a new chain.
#'
#' @param x An `nsbi_nre_posterior` object (named `x` to satisfy the [sample()]
#'   generic).
#' @param size,n Number of posterior draws (`n` is an alias for `size`).
#' @param obs Observation to condition on (defaults to the posterior's `x_obs`).
#' @param refresh Force a new run even when a cached one would do.
#' @param verbose Report sampling progress.
#' @param ... Unused.
#' @return An `n x dim` matrix of posterior draws (class `nsbi_samples`), with
#'   convergence diagnostics attached as attribute `diagnostics`.
#' @method sample nsbi_nre_posterior
#' @export
sample.nsbi_nre_posterior <- function(x, size = 1000, n = size, obs = NULL,
                                      refresh = FALSE, verbose = FALSE, ...) {
  mcmc_draws(x, n, obs, refresh, verbose)
}

#' @rdname log_prob
#' @export
log_prob.nsbi_nre_posterior <- function(post, theta, x = NULL,
                                        normalize = TRUE, ...) {
  mcmc_log_prob(post, theta, x, !missing(normalize) && isTRUE(normalize), "NRE")
}

#' @export
print.nsbi_nre_posterior <- function(x, ...) {
  cat_mcmc_posterior(x, "nsbi_nre_posterior")
  cat("  log_prob() is unnormalized: the evidence p(x) is not available.\n")
  cat("  sample(post, n), map_estimate(post)\n")
  invisible(x)
}
