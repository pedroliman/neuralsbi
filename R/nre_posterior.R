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
#'   from it, drawing more pools as needed if the first one lands fewer than
#'   `n_chains` draws in the posterior's support; `"proposal"` skips the
#'   weighting and keeps whichever prior draws land inside the posterior's
#'   support, also drawing more pools as needed. `"proposal"` is cheaper per
#'   accepted draw, but every draw the posterior excludes is wasted: if only a
#'   fraction `a` of the prior's mass survives, roughly `1 / a` draws are
#'   needed per starting point, so a posterior that rules out most of the
#'   prior can make `"proposal"` far slower than
#'   `"resample"` to find enough starting points, and it errors out once its
#'   draw budget is spent.
#' @param seed Optional integer seed.
#' @param max_batch Largest number of `(theta, x)` pairs evaluated in one call
#'   to the classifier, both at every MCMC step and in [log_prob()]. Only
#'   affects memory and speed; see `max_batch` on [log_ratio()] for the same
#'   knob on a direct ratio evaluation, and matters most here since the ratio
#'   classifier has no i.i.d. fast path -- every pair costs a forward pass.
#' @param ... Further arguments to the sampler: `width`, `max_steps` and
#'   `n_pool`.
#'
#' @return An object of class
#'   `c("nsbi_nre_posterior", "nsbi_mcmc_posterior", "nsbi_posterior")`.
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
                               seed = NULL, max_batch = 1e5, ...) {
  mcmc_posterior(fit, x_obs, "slice", n_chains, warmup, thin,
                 match.arg(init_strategy), seed, max_batch, list(...),
                 "nsbi_nre_posterior")
}
