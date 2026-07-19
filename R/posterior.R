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

#' Build a posterior from an NPE fit
#'
#' @param fit An `nsbi_npe` object from [npe()].
#' @param x_obs Optional default observation to condition on. If supplied it
#'   becomes the default `x` for [sample()], [log_prob()] and [map_estimate()].
#' @return An `nsbi_posterior` object.
#' @export
posterior <- function(fit, x_obs = NULL) {
  stopifnot(inherits(fit, "nsbi_npe"))
  if (!is.null(x_obs)) x_obs <- as_theta_matrix(x_obs, fit$dim_x)
  structure(
    list(fit = fit, x_obs = x_obs),
    class = "nsbi_posterior"
  )
}

#' @keywords internal
resolve_x <- function(post, x) {
  x <- x %||% post$x_obs
  if (is.null(x)) {
    stop("No observation supplied. Pass `x = ...` or set `x_obs` in posterior().",
         call. = FALSE)
  }
  as_theta_matrix(x, post$fit$dim_x)[1, , drop = FALSE]
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
  xo_std <- standardized_obs(post, obs)
  fit <- post$fit
  prior <- fit$prior
  bounded <- !is.null(prior$lower) || !is.null(prior$upper)

  collected <- matrix(0, nrow = 0, ncol = fit$dim_theta)
  n_tried <- 0L
  batch <- 0L
  while (nrow(collected) < n && batch < max_sampling_batches) {
    batch <- batch + 1L
    draw_std <- de_sample(fit$de, xo_std, n)
    draw <- invert_standardizer(fit$std_theta, draw_std)
    n_tried <- n_tried + n
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
#' @return Numeric vector of log posterior densities.
#' @export
log_prob <- function(post, theta, x = NULL, normalize = TRUE,
                     n_normalization = 10000L) {
  stopifnot(inherits(post, "nsbi_posterior"))
  fit <- post$fit
  theta <- as_theta_matrix(theta, fit$dim_theta)
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
#' @param post An `nsbi_posterior` object.
#' @param x Observation to condition on (defaults to `x_obs`).
#' @param n_init Number of initial draws used to seed the search.
#' @return Numeric vector: the MAP parameter estimate.
#' @export
map_estimate <- function(post, x = NULL, n_init = 1000L) {
  stopifnot(inherits(post, "nsbi_posterior"))
  fit <- post$fit
  draws <- sample.nsbi_posterior(post, n = n_init, obs = x)
  lp <- log_prob(post, draws, x = x, normalize = FALSE)
  start <- draws[which.max(lp), ]
  neg <- function(par) -log_prob(post, matrix(par, nrow = 1), x = x,
                                 normalize = FALSE)
  opt <- stats::optim(start, neg, method = "Nelder-Mead")
  opt$par
}

#' @export
print.nsbi_posterior <- function(x, ...) {
  cat("<nsbi_posterior>\n")
  cat(sprintf("  parameters (dim): %d\n", x$fit$dim_theta))
  cat(sprintf("  conditioned on x: %s\n",
              if (is.null(x$x_obs)) "(none set)" else
                paste(signif(x$x_obs[1, ], 4), collapse = ", ")))
  cat("  sample(post, n), log_prob(post, theta), map_estimate(post)\n")
  invisible(x)
}
