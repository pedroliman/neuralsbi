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
#' @param thin Keep one draw in `thin`. The default is 2: the slice width is
#'   adapted during warmup, and with that the retained draws are already close
#'   to independent -- on a Gaussian target with 20 chains, `thin = 2` gives a
#'   bulk ESS of about 96% of them. Since each evaluation is a forward pass over
#'   every observation, thinning harder buys very little for what it costs.
#'   Raise it if the reported ESS says you need to. (Python `sbi` thins by 1 by
#'   default; it thinned by 10 up to v0.21.)
#' @param init_strategy `"resample"` (default, and `sbi`'s) weights a pool of
#'   prior draws by the posterior density and resamples the starting points from
#'   it; `"proposal"` skips the weighting and keeps whichever prior draws land
#'   inside the posterior's support, drawing more pools as needed. `"proposal"`
#'   is cheaper per accepted draw, but every draw the posterior excludes is
#'   wasted: if only a fraction `a` of the prior's mass survives, roughly `1 /
#'   a` draws are needed per starting point, so a posterior that rules out most
#'   of the prior can make `"proposal"` far slower than `"resample"` to find
#'   enough starting points, and it errors out once its draw budget is spent.
#'   The pool is 1000 draws, set with `n_pool`. `sbi` uses 10,000 *per chain*,
#'   which on a surrogate summed over thousands of observations is a large bill
#'   before the first MCMC step.
#' @param seed Optional integer seed.
#' @param ... Further arguments to the sampler: `width`, `max_steps` and
#'   `n_pool` for `"slice"`, or `iter_warmup`, `iter_sampling` and `refresh`
#'   for `"stan"`.
#'
#' @return An object of class
#'   `c("nsbi_nle_posterior", "nsbi_mcmc_posterior", "nsbi_posterior")`.
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
  sampler <- match.arg(sampler)
  mcmc_posterior(fit, x_obs, sampler,
                 n_chains %||% if (sampler == "stan") 4L else 20L,
                 warmup, thin, match.arg(init_strategy), seed, list(...),
                 "nsbi_nle_posterior")
}

#' Check the arguments of an MCMC-sampled posterior and assemble it
#'
#' Shared by [posterior.nsbi_nle()] and [posterior.nsbi_nre()]. Both wrap a fit
#' and a sampler configuration around a draw cache, and every argument except
#' the sampler is checked the same way; only the class they carry and the
#' samplers they allow differ, and both of those are settled by the caller
#' before it gets here.
#'
#' @inheritParams posterior.nsbi_nle
#' @param dots The sampler arguments the caller collected from `...`.
#' @param class The posterior class to stamp on the result.
#' @keywords internal
mcmc_posterior <- function(fit, x_obs, sampler, n_chains, warmup, thin,
                           init_strategy, seed, dots, class) {
  check_fit_alive(fit)
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
  structure(
    list(
      fit = fit,
      x_obs = x_obs,
      sampler = sampler,
      control = list(n_chains = n_chains,
                     warmup = warmup,
                     thin = thin,
                     init_strategy = init_strategy,
                     seed = seed,
                     dots = dots),
      cache = new.env(parent = emptyenv())
    ),
    class = c(class, "nsbi_mcmc_posterior", "nsbi_posterior")
  )
}

#' Validate one MCMC count argument and return it as an integer
#'
#' `as.integer()` on its own lets a nonsensical count through. `thin = 0` is the
#' one that bites: [slice_sample()] then runs `warmup` iterations and no kept
#' ones, and returns its zero-initialized array of draws with diagnostics
#' computed on it. The counts are checked here, before anything is stored on the
#' posterior object.
#'
#' @param value The supplied value.
#' @param name Argument name, for the error message.
#' @param min Smallest value allowed.
#' @param why Optional clause explaining the bound.
#' @keywords internal
check_mcmc_count <- function(value, name, min, why = NULL) {
  ok <- is.numeric(value) && length(value) == 1L && is.finite(value) &&
    value == trunc(value) && value >= min
  if (!ok) {
    stop(sprintf("`%s` must be a single whole number, at least %d%s.",
                 name, min, if (is.null(why)) "" else paste0(", ", why)),
         call. = FALSE)
  }
  as.integer(value)
}

#' All rows of the observation, unlike `resolve_x()` which keeps only the first
#'
#' Thin wrapper around [resolve_obs()] with `first_row = FALSE`; see there for
#' the rationale shared with [resolve_x()], including why a non-finite entry
#' stops here rather than warns.
#'
#' @param post An `nsbi_mcmc_posterior` object.
#' @param x Observation to condition on, or `NULL` to use the posterior's
#'   `x_obs`.
#' @param arg Name the caller's argument goes by, for the error message.
#'   [sample()] calls it `obs` and [log_prob()] calls it `x`.
#' @keywords internal
resolve_x_iid <- function(post, x, arg = "obs") {
  resolve_obs(post, x, first_row = FALSE, arg = arg)
}

#' Sample an MCMC posterior
#'
#' Runs the sampler chosen by [posterior.nsbi_nle()] or [posterior.nsbi_nre()].
#' Draws are cached on the posterior object, so asking for the same or fewer
#' draws from the same observation returns immediately instead of re-running a
#' chain -- which is what makes [summary()] and repeated calls tolerable.
#'
#' @param x An `nsbi_mcmc_posterior` object (named `x` to satisfy the
#'   [sample()] generic).
#' @param size,n Number of posterior draws (`n` is an alias for `size`).
#' @param obs Observation to condition on (defaults to the posterior's `x_obs`).
#' @param refresh Force a new run even when a cached one would do.
#' @param verbose Report sampling progress.
#' @param ... Unused.
#' @return An `n x dim` matrix of posterior draws (class `nsbi_samples`), with
#'   convergence diagnostics attached as attribute `diagnostics`.
#' @method sample nsbi_mcmc_posterior
#' @export
sample.nsbi_mcmc_posterior <- function(x, size = 1000, n = size, obs = NULL,
                                       refresh = FALSE, verbose = FALSE, ...) {
  n <- check_count(n, "n", why = "since it is the number of posterior draws")
  mcmc_draws(x, n, obs, refresh, verbose)
}

#' Run (or reuse) the chain behind an MCMC posterior's [sample()] method
#'
#' The body [sample.nsbi_mcmc_posterior()] runs regardless of whether the
#' underlying fit is an [nle()] or [nre()]. Which sampler runs is read off the
#' posterior object, so nothing here needs to know which kind of fit produced
#' it.
#'
#' @param post An MCMC-sampled `nsbi_posterior`.
#' @param n Number of draws.
#' @param obs Observation to condition on, or `NULL` for the posterior's own.
#' @param refresh Force a new run even when a cached one would do.
#' @param verbose Report sampling progress.
#' @keywords internal
mcmc_draws <- function(post, n, obs, refresh, verbose) {
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
    slice_sample_surrogate(fit, x_obs, ctl, n, verbose = verbose)
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

#' Slice-sample the unnormalized posterior of an [nle()] or [nre()] fit
#'
#' The sampler does not care which surrogate produced the potential, so
#' [surrogate_potential()] is the only line that looks at which fit it has.
#' @keywords internal
slice_sample_surrogate <- function(fit, x_obs, ctl, n, verbose = FALSE) {
  dots <- ctl$dots
  potential <- surrogate_potential(fit, x_obs)
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
  diagnostics <- mcmc_diagnostics(res$chains)
  # n_evals is the cost of the run -- it climbs fastest exactly when the slice
  # width is mismatched to the target, which is what warmup adapts it for (see
  # ?nsbi_mcmc). Carried as an attribute rather than a column since it is one
  # number for the whole run, not one per parameter like rhat and ess_bulk.
  attr(diagnostics, "n_evals") <- res$n_evals
  list(draws = res$draws, diagnostics = diagnostics)
}

#' A per-coordinate scale for the prior, used to size the initial slice width
#' @keywords internal
prior_scale <- function(prior) {
  if (!is.null(prior$lower) && !is.null(prior$upper) &&
      all(is.finite(prior$lower)) && all(is.finite(prior$upper))) {
    return(pmax((prior$upper - prior$lower) / 10, .Machine$double.eps))
  }
  s <- apply(sample_prior(prior, 1000L), 2, stats::sd)
  s[!is.finite(s) | s <= 0] <- 1
  s
}

#' @rdname log_prob
#' @export
log_prob.nsbi_mcmc_posterior <- function(post, theta, x = NULL,
                                         normalize = TRUE, ...) {
  what <- if (inherits(post$fit, "nsbi_nle")) "NLE" else "NRE"
  mcmc_log_prob(post, theta, x, !missing(normalize) && isTRUE(normalize), what)
}

#' The unnormalized log density behind an MCMC posterior's [log_prob()] method
#'
#' Neither surrogate gives the evidence \eqn{p(x)}, so `normalize = TRUE` has
#' nothing to normalize by and says so rather than returning a number that
#' looks like a density.
#'
#' @param post An MCMC-sampled `nsbi_posterior`.
#' @param theta Parameter values to evaluate.
#' @param x Observation to condition on, or `NULL` for the posterior's own.
#' @param warn Warn that `normalize` is being ignored (the caller decides,
#'   since only it can see whether the argument was actually supplied).
#' @param what The method name to use in that warning, `"NLE"` or `"NRE"`.
#' @keywords internal
mcmc_log_prob <- function(post, theta, x, warn, what) {
  if (isTRUE(warn)) {
    warning(sprintf(paste0("An %s posterior has no normalizing constant; ",
                           "`normalize` is ignored and the value returned is ",
                           "unnormalized."), what), call. = FALSE)
  }
  # surrogate_potential()'s closure only reshapes theta with
  # as_theta_matrix(), which coerces a non-numeric column to NA rather than
  # erroring -- so a bad theta used to come back as a silent NA log-prob
  # instead of the same named error surrogate_score() already gives
  # log_lik()/log_ratio() for the same fit (#163). Checked here, once, with
  # the caller's theta in hand, rather than inside surrogate_potential()'s
  # returned closure, which every MCMC step also calls with a theta it built
  # internally and has no need to re-validate.
  theta <- check_matrix(theta, post$fit$dim_theta, "theta",
                        "one parameter per column")
  potential <- surrogate_potential(post$fit, resolve_x_iid(post, x, "x"))
  potential(theta)
}

#' @export
print.nsbi_mcmc_posterior <- function(x, ...) {
  cat_mcmc_posterior(x, class(x)[1])
  cat("  log_prob() is unnormalized: the evidence p(x) is not available.\n")
  if (inherits(x$fit, "nsbi_nle")) {
    cat("  sample(post, n), map_estimate(post), stan_code(post$fit)\n")
  } else {
    cat("  sample(post, n), map_estimate(post)\n")
  }
  invisible(x)
}

#' The summary block every MCMC posterior's `print()` method shares
#'
#' Everything except the closing "here is what you can do with it" lines, which
#' differ: an [nle()] fit can be exported to Stan and an [nre()] fit cannot.
#' @param x The posterior object.
#' @param class Class name to print in the header.
#' @keywords internal
cat_mcmc_posterior <- function(x, class) {
  cat(sprintf("<%s>\n", class))
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
    cat(sprintf("  last run        : %s\n",
                format_mcmc_diagnostics(x$cache$diagnostics)))
  }
  invisible(x)
}
