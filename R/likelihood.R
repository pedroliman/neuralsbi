#' Evaluate a surrogate likelihood
#'
#' `log_lik()` evaluates the likelihood learned by [nle()]:
#' \eqn{\log q_\phi(x \mid \theta)}, in the original units of both.
#'
#' Rows of `x` are treated as **independent observations from the same
#' parameter**, so by default the result sums over them,
#' \eqn{\sum_i \log q_\phi(x_i \mid \theta)}. That sum is the whole point of
#' NLE: the estimator is trained on one observation at a time, and the number
#' of observations you condition on afterwards is free.
#'
#' @param fit An `nsbi_nle` fit from [nle()].
#' @param theta Parameter values: a numeric vector (one parameter set) or an
#'   `n_theta x dim_theta` matrix.
#' @param x Observed data: a numeric vector (one observation) or an
#'   `n_obs x dim_x` matrix whose rows are independent observations.
#' @param sum_iid Sum the log-density over the rows of `x` (the default). Set
#'   `FALSE` to get the per-observation values instead.
#' @param max_batch Largest number of `(theta, x)` pairs evaluated in one call
#'   to the estimator. Only affects memory and speed.
#' @param ... Unused, for S3 consistency.
#'
#' @return With `sum_iid = TRUE`, a numeric vector with one entry per row of
#'   `theta`. With `sum_iid = FALSE`, an `n_theta x n_obs` matrix.
#'
#' @seealso [likelihood_fn()] for a closure over a fixed observation,
#'   [posterior()] to turn the likelihood into posterior draws.
#'
#' @examples
#' prior <- prior_uniform(c(mu = -3), c(mu = 3))
#' fit <- nle(prior, function(mu) c(y = rnorm(1, mu, 0.5)),
#'            n_simulations = 1000, density_estimator = "linear_gaussian")
#'
#' x_obs <- matrix(rnorm(20, mean = 1, sd = 0.5), ncol = 1)
#' grid <- matrix(seq(-2, 2, length.out = 5), ncol = 1)
#' log_lik(fit, grid, x_obs)
#' @export
log_lik <- function(fit, theta, x, ...) UseMethod("log_lik")

#' @rdname log_lik
#' @export
log_lik.nsbi_nle <- function(fit, theta, x, sum_iid = TRUE,
                             max_batch = 1e5, ...) {
  check_fit_alive(fit)
  # Both arguments are matrices with a required width, so a bare "expected 2
  # columns" would leave the user guessing which one is wrong.
  theta <- as_lik_matrix(theta, fit$dim_theta, "theta", "one parameter per column")
  x <- as_lik_matrix(x, fit$dim_x, "x",
                     "one row per independent observation")

  theta_z <- apply_standardizer(fit$std_theta, theta)
  x_z <- apply_standardizer(fit$std_x, x)
  n_theta <- nrow(theta_z)
  n_obs <- nrow(x_z)

  # The estimator broadcasts one conditioning row over many targets but not the
  # reverse, so the (theta, x) cross product is expanded explicitly. Blocking
  # over theta keeps the expansion bounded no matter how many trials there are.
  per_block <- max(1L, floor(max_batch / max(n_obs, 1L)))
  out <- matrix(0, nrow = n_theta, ncol = n_obs)
  starts <- seq.int(1L, n_theta, by = per_block)
  for (s in starts) {
    idx <- seq.int(s, min(s + per_block - 1L, n_theta))
    cond <- theta_z[rep(idx, each = n_obs), , drop = FALSE]
    target <- x_z[rep(seq_len(n_obs), times = length(idx)), , drop = FALSE]
    lp <- de_log_prob(fit$de, target, cond)
    out[idx, ] <- matrix(lp, nrow = length(idx), ncol = n_obs, byrow = TRUE)
  }

  # de_log_prob works in standardized x space; this puts it back in the units
  # the simulator returned.
  out <- out + standardizer_log_jac(fit$std_x)

  if (!isTRUE(sum_iid)) {
    return(out)
  }
  stats::setNames(rowSums(out), rownames(theta))
}

#' Coerce one of `log_lik()`'s two matrix arguments, naming it if it is wrong
#' @keywords internal
as_lik_matrix <- function(value, d, arg, what) {
  tryCatch(
    as_theta_matrix(value, d),
    error = function(e) {
      stop(sprintf("`%s` must have %d columns (%s).", arg, d, what),
           call. = FALSE)
    }
  )
}

#' A surrogate likelihood as a plain R function
#'
#' Fixes an observation and returns `function(theta)`, so the learned
#' likelihood can be handed to anything in R that wants a log-density:
#' [stats::optim()], an MCMC package, an importance sampler, a profile
#' likelihood. The returned function is vectorized -- pass a matrix of
#' parameter values and get one log-likelihood per row.
#'
#' @param fit An `nsbi_nle` fit from [nle()].
#' @param x_obs The observation to condition on. Rows are independent
#'   observations; the log-likelihood sums over them.
#' @param ... Passed to [log_lik()].
#'
#' @return A function of `theta` returning a numeric vector of
#'   log-likelihoods, with the observation attached as attribute `x_obs`.
#'
#' @examples
#' prior <- prior_uniform(c(mu = -3), c(mu = 3))
#' fit <- nle(prior, function(mu) c(y = rnorm(1, mu, 0.5)),
#'            n_simulations = 1000, density_estimator = "linear_gaussian")
#'
#' loglik <- likelihood_fn(fit, matrix(rnorm(20, 1, 0.5), ncol = 1))
#' optimize(function(m) loglik(m), c(-3, 3), maximum = TRUE)$maximum
#' @export
likelihood_fn <- function(fit, x_obs, ...) {
  stopifnot(inherits(fit, "nsbi_nle"))
  x_obs <- as_theta_matrix(x_obs, fit$dim_x)
  force(x_obs)
  dots <- list(...)
  f <- function(theta) {
    do.call(log_lik, c(list(fit = fit, theta = theta, x = x_obs), dots))
  }
  attr(f, "x_obs") <- x_obs
  f
}

#' Unnormalized log posterior of an NLE fit
#'
#' \eqn{\log q_\phi(x \mid \theta) + \log p(\theta)}, returning `-Inf` outside
#' the prior support. This is the potential the MCMC samplers target.
#' @keywords internal
nle_potential <- function(fit, x_obs, max_batch = 1e5) {
  prior <- fit$prior
  # prior_custom() without a log_prob_fn returns NA rather than nothing, so the
  # only way to find out is to ask it. Better here than as a puzzling
  # initialization failure a few hundred lines later.
  probe <- tryCatch(prior$log_prob(sample_prior(prior, 2L)),
                    error = function(e) NA_real_)
  if (is.null(prior$log_prob) || all(is.na(probe))) {
    stop("Sampling an NLE posterior needs a prior log-density, and this prior ",
         "does not have one.\nRebuild it with prior_custom(..., log_prob_fn = ).",
         call. = FALSE)
  }
  x_obs <- as_theta_matrix(x_obs, fit$dim_x)
  bounded <- !is.null(prior$lower) || !is.null(prior$upper)
  function(theta) {
    theta <- as_theta_matrix(theta, fit$dim_theta)
    lp <- as.numeric(prior$log_prob(theta))
    if (bounded) lp[!within_support(prior, theta)] <- -Inf
    ok <- is.finite(lp)
    if (!any(ok)) return(lp)
    # Only evaluate the network where the prior gives the point any mass; that
    # is a large saving once a chain starts probing a bounded edge.
    lp[ok] <- lp[ok] + log_lik(fit, theta[ok, , drop = FALSE], x_obs,
                               max_batch = max_batch)
    lp
  }
}
