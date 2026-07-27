#' Posterior from a neural likelihood
#'
#' Bayes' rule turns an [nle()] fit into a posterior,
#' \eqn{p(\theta \mid x) \propto q_\phi(x \mid \theta)\,p(\theta)}, but the
#' result has no closed form and no direct sampler. `posterior()` on an
#' `nsbi_nle` therefore returns an object that samples with MCMC, and
#' [sample()] on it runs a chain rather than a forward pass.
#'
#' Two samplers are available. `"slice"` is the default: a vectorized
#' univariate slice sampler (see [nsbi_mcmc]) with nothing to tune and no
#' dependencies. `"stan"` writes the fitted likelihood out as a Stan program
#' (see [stan_code()]) and runs NUTS on it through \pkg{cmdstanr} or
#' \pkg{rstan}, which mixes better on correlated posteriors at the cost of a
#' one-time model compile.
#'
#' @param fit An `nsbi_nle` object from [nle()].
#' @param x_obs Observation to condition on. Rows are treated as independent
#'   observations from the same parameter, and the log-likelihood sums over
#'   them.
#' @param sampler `"slice"` (the default) or `"stan"`.
#' @param n_chains Number of chains. The default depends on the sampler: 20 for
#'   `"slice"`, which evaluates every chain in one batched call per step and so
#'   pays almost nothing for more of them, and 4 for `"stan"`, where each chain
#'   is a separate process with its own warmup to pay for.
#' @param warmup Steps discarded at the start of each chain.
#' @param thin Keep one draw in `thin`. Python `sbi` defaults to 10; the
#'   default here is 2, because the slice width is adapted during warmup and
#'   the resulting draws are already close to independent -- on a Gaussian
#'   target with 20 chains, `thin = 2` gives a bulk ESS of about 96% of the
#'   draws. Since each evaluation is a forward pass over every observation,
#'   thinning by 10 would cost five times as much for that last 4%. Raise it if
#'   the reported ESS says you need to.
#' @param init_strategy `"resample"` (default) weights a pool of prior draws by
#'   the posterior density and resamples the starting points from it;
#'   `"proposal"` starts from plain prior draws.
#' @param seed Optional integer seed.
#' @param ... Further arguments to the sampler: `width`, `max_steps` and
#'   `n_pool` for `"slice"`, or `iter_warmup`, `iter_sampling` and `refresh`
#'   for `"stan"`.
#'
#' @return An object of class `c("nsbi_nle_posterior", "nsbi_posterior")`.
#'
#' @examples
#' prior <- prior_uniform(c(mu = -3), c(mu = 3))
#' fit <- nle(prior, function(mu) c(y = rnorm(1, mu, 0.5)),
#'            n_simulations = 1000, density_estimator = "linear_gaussian")
#'
#' x_obs <- matrix(rnorm(50, mean = 1, sd = 0.5), ncol = 1)
#' post <- posterior(fit, x_obs, n_chains = 4, warmup = 50, thin = 2)
#' draws <- sample(post, 400)
#' colMeans(draws)
#' @export
posterior.nsbi_nle <- function(fit, x_obs = NULL,
                               sampler = c("slice", "stan"),
                               n_chains = NULL, warmup = 200L, thin = 2L,
                               init_strategy = c("resample", "proposal"),
                               seed = NULL, ...) {
  check_fit_alive(fit)
  sampler <- match.arg(sampler)
  init_strategy <- match.arg(init_strategy)
  if (!is.null(x_obs)) x_obs <- as_theta_matrix(x_obs, fit$dim_x)
  n_chains <- n_chains %||% if (sampler == "stan") 4L else 20L
  if (n_chains < 2L) {
    stop("`n_chains` must be at least 2, so convergence can be diagnosed.",
         call. = FALSE)
  }

  structure(
    list(
      fit = fit,
      x_obs = x_obs,
      sampler = sampler,
      control = list(n_chains = as.integer(n_chains),
                     warmup = as.integer(warmup),
                     thin = as.integer(thin),
                     init_strategy = init_strategy,
                     seed = seed,
                     dots = list(...)),
      cache = new.env(parent = emptyenv())
    ),
    class = c("nsbi_nle_posterior", "nsbi_posterior")
  )
}

#' All rows of the observation, unlike [resolve_x()] which keeps only the first
#' @keywords internal
resolve_x_iid <- function(post, x) {
  x <- x %||% post$x_obs
  if (is.null(x)) {
    stop("No observation supplied. Pass `obs = ...` or set `x_obs` in posterior().",
         call. = FALSE)
  }
  as_theta_matrix(x, post$fit$dim_x)
}

#' Sample an NLE posterior with MCMC
#'
#' Runs the sampler chosen by [posterior.nsbi_nle()]. Draws are cached on the
#' posterior object, so asking for the same or fewer draws from the same
#' observation returns immediately instead of re-running a chain -- which is
#' what makes [summary()] and repeated calls tolerable.
#'
#' @param x An `nsbi_nle_posterior` object (named `x` to satisfy the [sample()]
#'   generic).
#' @param size,n Number of posterior draws (`n` is an alias for `size`).
#' @param obs Observation to condition on (defaults to the posterior's `x_obs`).
#' @param refresh Force a new run even when a cached one would do.
#' @param verbose Report sampling progress.
#' @param ... Unused.
#' @return An `n x dim` matrix of posterior draws (class `nsbi_samples`), with
#'   convergence diagnostics attached as attribute `diagnostics`.
#' @method sample nsbi_nle_posterior
#' @export
sample.nsbi_nle_posterior <- function(x, size = 1000, n = size, obs = NULL,
                                      refresh = FALSE, verbose = FALSE, ...) {
  post <- x
  fit <- post$fit
  x_obs <- resolve_x_iid(post, obs)

  cached <- post$cache$draws
  if (!isTRUE(refresh) && !is.null(cached) &&
      identical(post$cache$x_obs, x_obs) && nrow(cached) >= n) {
    return(finish_draws(cached[seq_len(n), , drop = FALSE],
                        post$cache$diagnostics, fit))
  }

  ctl <- post$control
  if (!is.null(ctl$seed)) set.seed(ctl$seed)

  run <- if (post$sampler == "stan") {
    stan_sample_nle(fit, x_obs, ctl, n, verbose = verbose)
  } else {
    slice_sample_nle(fit, x_obs, ctl, n, verbose = verbose)
  }

  post$cache$draws <- run$draws
  post$cache$diagnostics <- run$diagnostics
  post$cache$x_obs <- x_obs
  finish_draws(run$draws[seq_len(min(n, nrow(run$draws))), , drop = FALSE],
               run$diagnostics, fit)
}

#' @keywords internal
finish_draws <- function(draws, diagnostics, fit) {
  if (is.null(colnames(draws)) && !is.null(fit$param_names)) {
    colnames(draws) <- fit$param_names
  }
  if (!is.null(diagnostics) && !is.null(fit$param_names)) {
    rownames(diagnostics) <- fit$param_names
  }
  attr(draws, "diagnostics") <- diagnostics
  structure(draws, class = c("nsbi_samples", class(draws)))
}

#' @keywords internal
slice_sample_nle <- function(fit, x_obs, ctl, n, verbose = FALSE) {
  dots <- ctl$dots
  potential <- nle_potential(fit, x_obs)
  init <- mcmc_init(fit$prior, potential, ctl$n_chains,
                    strategy = ctl$init_strategy,
                    n_pool = dots$n_pool %||% 1000L)

  # Scale the initial slice width to the prior, so the sampler starts with the
  # right order of magnitude in each coordinate rather than a bare 1.
  width <- dots$width %||% prior_scale(fit$prior)

  res <- slice_sample(potential, init, n_draws = n,
                      warmup = ctl$warmup, thin = ctl$thin,
                      width = width, max_steps = dots$max_steps %||% 100L,
                      verbose = verbose)
  list(draws = res$draws, diagnostics = mcmc_diagnostics(res$chains))
}

#' A per-coordinate scale for the prior, used to size the initial slice width
#' @keywords internal
prior_scale <- function(prior, n = 1000L) {
  if (!is.null(prior$lower) && !is.null(prior$upper) &&
      all(is.finite(prior$lower)) && all(is.finite(prior$upper))) {
    return(pmax((prior$upper - prior$lower) / 10, .Machine$double.eps))
  }
  s <- apply(sample_prior(prior, n), 2, stats::sd)
  s[!is.finite(s) | s <= 0] <- 1
  s
}

#' @rdname log_prob
#' @export
log_prob.nsbi_nle_posterior <- function(post, theta, x = NULL,
                                        normalize = TRUE, ...) {
  if (!missing(normalize) && isTRUE(normalize)) {
    warning("An NLE posterior has no normalizing constant; `normalize` is ",
            "ignored and the value returned is unnormalized.", call. = FALSE)
  }
  potential <- nle_potential(post$fit, resolve_x_iid(post, x))
  potential(theta)
}

#' @export
print.nsbi_nle_posterior <- function(x, ...) {
  cat("<nsbi_nle_posterior>\n")
  cat(sprintf("  parameters (dim): %d\n", x$fit$dim_theta))
  if (!is.null(x$fit$param_names)) {
    cat("    names         :", paste(x$fit$param_names, collapse = ", "), "\n")
  }
  cat(sprintf("  observations    : %s\n",
              if (is.null(x$x_obs)) "(none set)"
              else sprintf("%d x %d", nrow(x$x_obs), ncol(x$x_obs))))
  ctl <- x$control
  cat(sprintf("  sampler         : %s (%d chains, warmup %d, thin %d, init %s)\n",
              x$sampler, ctl$n_chains, ctl$warmup, ctl$thin, ctl$init_strategy))
  if (!is.null(x$cache$diagnostics)) {
    d <- x$cache$diagnostics
    cat(sprintf("  last run        : max Rhat %.3f, min bulk ESS %.0f\n",
                max(d$rhat, na.rm = TRUE), min(d$ess_bulk, na.rm = TRUE)))
  }
  cat("  log_prob() is unnormalized: the evidence p(x) is not available.\n")
  cat("  sample(post, n), map_estimate(post), stan_code(post$fit)\n")
  invisible(x)
}
