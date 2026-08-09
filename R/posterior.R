#' Posterior objects
#'
#' A posterior wraps a trained [npe()] fit together with (optionally) a default
#' observation `x_obs`. It knows how to draw posterior samples, evaluate the
#' posterior log-density, and find the maximum-a-posteriori (MAP) estimate. All
#' transforms between standardized training space and the original parameter
#' space are handled internally.
#'
#' For bounded priors, samples that fall outside the prior support are rejected
#' ("leakage" correction), and [log_prob()] is renormalized by the estimated
#' acceptance probability so it integrates to one over the support.
#'
#' @name posterior
NULL

#' Build a posterior from a fit
#'
#' The two inference methods reach a posterior by different routes, and
#' `posterior()` hides the difference. An [npe()] fit already *is* a posterior
#' estimator, so the returned object samples with a forward pass. An [nle()]
#' fit only knows the likelihood, so the returned object samples with MCMC and
#' takes the extra arguments that implies.
#'
#' @param fit An `nsbi_npe` object from [npe()] or [npe_sequential()], or an
#'   `nsbi_nle` object from [nle()].
#' @param x_obs Optional default observation to condition on. If supplied it
#'   becomes the default `x` for [sample()], [log_prob()] and [map_estimate()].
#'   For an NLE fit, rows of `x_obs` are independent observations.
#' @param ... Passed to methods. See [posterior.nsbi_nle()] for the MCMC
#'   controls an NLE fit accepts.
#' @return An `nsbi_posterior` object.
#' @seealso [save_npe()], which is how a torch-backed fit gets to disk and
#'   back; `readRDS()` returns one whose network is dead, and `posterior()`
#'   says so rather than failing later.
#' @export
posterior <- function(fit, x_obs = NULL, ...) UseMethod("posterior")

#' @rdname posterior
#' @export
posterior.default <- function(fit, x_obs = NULL, ...) {
  stop("posterior() needs a fit from npe() or nle(), not an object of class ",
       paste(class(fit), collapse = "/"), ".", call. = FALSE)
}

#' @rdname posterior
#' @export
posterior.nsbi_npe <- function(fit, x_obs = NULL, ...) {
  check_fit_alive(fit)
  if (!is.null(x_obs)) {
    x_obs <- check_numeric(x_obs, "x_obs")
    check_finite(x_obs, "x_obs")
    x_obs <- as_theta_matrix(x_obs, fit$dim_x)
  }
  structure(
    list(fit = fit, x_obs = x_obs),
    class = "nsbi_posterior"
  )
}

#' Resolve the observation a posterior conditions on
#'
#' `resolve_x()` and `resolve_x_iid()` are both thin wrappers around this: an
#' NPE fit maps one observation to one posterior, so `resolve_x()` calls this
#' with `first_row = TRUE` and keeps only row 1; an NLE fit's log-likelihood
#' sums over rows as independent observations of the same parameter, so
#' `resolve_x_iid()` calls this with `first_row = FALSE` and keeps every row.
#' The same `x_obs` therefore means "the first observation" to [npe()] and
#' "200 observations" to [nle()], and a user moving a working call from one to
#' the other would otherwise get a posterior conditioned on a single data
#' point with nothing said about it. Truncating to row 1 warns rather than
#' fails: taking row 1 of a simulation matrix is a reasonable thing to ask
#' for.
#'
#' A non-finite entry is a different matter and stops here regardless of
#' `first_row`. An `NA` passes through `apply_standardizer()` into
#' `de_sample()` (NPE) or into [mcmc_init()] (NLE) and comes back as
#' unusable output with nothing to say about the observation -- for NPE the
#' complaint lands in `stats::quantile()` inside `summary()`, for NLE in
#' `mcmc_init()` complaining about initialization. There is nothing sensible
#' to condition on either way, so this errors rather than warns.
#'
#' @param post An `nsbi_posterior` object.
#' @param x Observation to condition on, or `NULL` to use the posterior's
#'   `x_obs`.
#' @param first_row `TRUE` to keep row 1 only, with a warning when there was
#'   more than one row (the [resolve_x()] behavior); `FALSE` to keep every row
#'   (the [resolve_x_iid()] behavior).
#' @param arg Name the caller's argument goes by, for `check_numeric()`'s and
#'   `check_finite()`'s error messages. Defaults to `"x"` when `first_row` and
#'   `"obs"` otherwise, matching [sample()]'s and [log_prob()]'s parameter
#'   names; the "no observation supplied" message is fixed to the same names
#'   and does not vary with `arg`.
#' @return A matrix with `dim_x` columns: one row if `first_row`, every row
#'   otherwise.
#' @keywords internal
resolve_obs <- function(post, x, first_row, arg = if (first_row) "x" else "obs") {
  check_arg <- if (is.null(x)) "x_obs" else arg
  x <- x %||% post$x_obs
  if (is.null(x)) {
    stop(sprintf(
      "No observation supplied. Pass `%s = ...` or set `x_obs` in posterior().",
      if (first_row) "x" else "obs"),
      call. = FALSE)
  }
  x <- check_numeric(x, check_arg)
  check_finite(x, check_arg)
  x <- as_theta_matrix(x, post$fit$dim_x)
  if (!first_row) return(x)
  if (nrow(x) > 1L) {
    warning(sprintf(
      "Observation has %d rows; an NPE posterior conditions on one, so only row 1 is used. ",
      nrow(x)),
      "Use nle() if the rows are repeated observations of the same parameter.",
      call. = FALSE)
  }
  x[1, , drop = FALSE]
}

#' The single observation an NPE posterior conditions on
#'
#' Thin wrapper around [resolve_obs()] with `first_row = TRUE`; see there for
#' the rationale shared with [resolve_x_iid()].
#'
#' @param post An `nsbi_posterior` object.
#' @param x Observation to condition on, or `NULL` to use the posterior's
#'   `x_obs`.
#' @return A one-row matrix with `dim_x` columns.
#' @keywords internal
resolve_x <- function(post, x) {
  resolve_obs(post, x, first_row = TRUE)
}

#' Sample from a posterior
#'
#' @param x An `nsbi_posterior` object (named `x` to satisfy the [sample()]
#'   generic).
#' @param size,n Number of posterior draws (`n` is an alias for `size`).
#' @param obs Observation to condition on (defaults to the posterior's `x_obs`).
#' @param max_sampling_batches Safety cap on rejection-sampling rounds for
#'   bounded priors.
#' @param ... Unused.
#' @return An `n x dim` matrix of posterior draws (class `nsbi_samples`).
#' @method sample nsbi_posterior
#' @export
sample.nsbi_posterior <- function(x, size = 1000, n = size, obs = NULL,
                                  max_sampling_batches = 100L, ...) {
  post <- x
  max_sampling_batches <- check_count(
    max_sampling_batches, "max_sampling_batches",
    why = "since one batch is one round of rejection sampling")
  xo_std <- standardized_obs(post, obs)
  fit <- post$fit
  prior <- fit$prior
  bounded <- !is.null(prior$lower) || !is.null(prior$upper)

  collected <- matrix(0, nrow = 0, ncol = fit$dim_theta)
  n_tried <- 0L
  batch <- 0L
  while (nrow(collected) < n && batch < max_sampling_batches) {
    batch <- batch + 1L
    n_needed <- n - nrow(collected)
    draw_std <- de_sample(fit$de, xo_std, n_needed)
    draw <- invert_standardizer(fit$std_theta, draw_std)
    n_tried <- n_tried + n_needed
    if (bounded) {
      draw <- draw[within_support(prior, draw), , drop = FALSE]
    }
    collected <- rbind(collected, draw)
  }
  if (nrow(collected) < n) {
    warning(sprintf(
      "Only %d/%d samples inside prior support after %d batches (acceptance %.2f). ",
      nrow(collected), n, batch, nrow(collected) / max(n_tried, 1)),
      "The estimator is leaking mass outside the prior; consider more simulations.",
      call. = FALSE)
  }
  out <- collected[seq_len(min(n, nrow(collected))), , drop = FALSE]
  if (is.null(colnames(out)) && !is.null(fit$param_names)) {
    colnames(out) <- fit$param_names
  }
  attr(out, "acceptance_rate") <- nrow(collected) / max(n_tried, 1)
  structure(out, class = c("nsbi_samples", class(out)))
}

#' @keywords internal
standardized_obs <- function(post, obs) {
  xo <- resolve_x(post, obs)
  apply_standardizer(post$fit$std_x, xo)
}

#' Posterior log-density
#'
#' @param post An `nsbi_posterior` object.
#' @param theta Matrix (or vector) of parameter values to evaluate.
#' @param x Observation to condition on (defaults to `x_obs`).
#' @param normalize For bounded priors, renormalize by the estimated acceptance
#'   probability and return `-Inf` outside the prior support.
#' @param n_normalization Number of draws used to estimate the normalizing
#'   (acceptance) constant when `normalize = TRUE`.
#' @param ... Passed to methods.
#' @return Numeric vector of log posterior densities. For a posterior built
#'   from an [nle()] fit the value is **unnormalized** -- the evidence
#'   \eqn{p(x)} is not available -- so differences between two `theta` are
#'   meaningful but the absolute level is not.
#' @export
log_prob <- function(post, theta, x = NULL, ...) UseMethod("log_prob")

#' @rdname log_prob
#' @export
log_prob.nsbi_posterior <- function(post, theta, x = NULL, normalize = TRUE,
                                    n_normalization = 10000L, ...) {
  n_normalization <- check_count(n_normalization, "n_normalization")
  fit <- post$fit
  theta <- as_theta_matrix(check_numeric(theta, "theta"), fit$dim_theta)
  xo_std <- standardized_obs(post, x)
  theta_z <- apply_standardizer(fit$std_theta, theta)
  lp <- de_log_prob(fit$de, theta_z, xo_std) + standardizer_log_jac(fit$std_theta)

  prior <- fit$prior
  bounded <- !is.null(prior$lower) || !is.null(prior$upper)
  if (normalize && bounded) {
    draw_std <- de_sample(fit$de, xo_std, n_normalization)
    draw <- invert_standardizer(fit$std_theta, draw_std)
    acc <- mean(within_support(prior, draw))
    acc <- max(acc, 1 / n_normalization)
    lp <- lp - log(acc)
    lp[!within_support(prior, theta)] <- -Inf
  }
  lp
}

#' Maximum a posteriori (MAP) estimate
#'
#' Starts from the best of a set of posterior draws and refines with a
#' derivative-free optimizer.
#'
#' On a posterior from [nle()] the initial draws come from MCMC, so `n_init`
#' buys a chain rather than a forward pass. They are cached on the posterior
#' like any other run.
#'
#' @param post An `nsbi_posterior` object.
#' @param x Observation to condition on (defaults to `x_obs`).
#' @param n_init Number of initial draws used to seed the search.
#' @return Numeric vector: the MAP parameter estimate.
#' @export
map_estimate <- function(post, x = NULL, n_init = 1000L) {
  stopifnot(inherits(post, "nsbi_posterior"))
  n_init <- check_count(n_init, "n_init",
                        why = "since the search starts from the best of them")
  fit <- post$fit
  # Through the generic, not sample.nsbi_posterior() directly: an nle() fit's
  # estimator has the roles swapped, so the NPE sampler asked to draw from it
  # returns draws in x space. Where dim_x and dim_theta differ that is an error
  # about non-conformable arguments; where they happen to match it is a set of
  # starting points quietly drawn from the wrong distribution.
  draws <- sample(post, n = n_init, obs = x)
  lp <- log_prob(post, draws, x = x, normalize = FALSE)
  start <- draws[which.max(lp), ]
  neg <- function(par) -log_prob(post, matrix(par, nrow = 1), x = x,
                                 normalize = FALSE)
  opt <- stats::optim(start, neg, method = "Nelder-Mead")
  out <- opt$par
  if (is.null(names(out))) names(out) <- fit$param_names
  out
}

#' @export
print.nsbi_posterior <- function(x, ...) {
  cat("<nsbi_posterior>\n")
  cat(sprintf("  parameters (dim): %d\n", x$fit$dim_theta))
  if (!is.null(x$fit$param_names)) {
    cat("    names         :", paste(x$fit$param_names, collapse = ", "), "\n")
  }
  cat(sprintf("  conditioned on x: %s\n",
              if (is.null(x$x_obs)) "(none set)" else {
                vals <- signif(x$x_obs[1, ], 4)
                nm <- x$fit$x_names
                if (!is.null(nm)) vals <- stats::setNames(vals, nm)
                paste(if (!is.null(nm)) paste0(nm, "=", vals) else vals,
                     collapse = ", ")
              }))
  cat("  sample(post, n), log_prob(post, theta), map_estimate(post)\n")
  invisible(x)
}
