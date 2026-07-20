#' Build a posterior from a trained NPE fit
#'
#' Wraps a trained [NPE] object as a [DirectPosterior]: the density estimator
#' *is* the posterior, so sampling and density evaluation are immediate. All
#' standardization transforms are handled internally; for bounded priors,
#' out-of-support draws are rejected and [log_prob()] is renormalized by the
#' estimated acceptance probability.
#'
#' @param fit A trained [NPE] object from [npe()].
#' @param x_obs Optional default observation to condition on.
#' @return A [DirectPosterior] object.
#' @export
posterior <- function(fit, x_obs = NULL) {
  if (!S7_inherits(fit, NPE)) stop("`fit` must be an NPE object.", call. = FALSE)
  build_posterior(fit, x_obs = x_obs)
}

#' @keywords internal
resolve_x <- function(post, x) {
  x <- x %||% post@default_x
  if (is.null(x)) {
    stop("No observation supplied. Pass `x = ...` or set `x_obs` in posterior().",
         call. = FALSE)
  }
  as_theta_matrix(x, ncol(post@fit@x))[1, , drop = FALSE]
}

#' @keywords internal
standardized_obs <- function(post, obs) {
  apply_standardizer(post@fit@std_x, resolve_x(post, obs))
}

method(sample, DirectPosterior) <- function(x, size = 1000, n = size,
                                            obs = NULL,
                                            max_sampling_batches = 100L, ...) {
  post <- x
  fit <- post@fit
  xo_std <- standardized_obs(post, obs)
  prior <- fit@prior
  bounded <- !is.null(prior@lower) || !is.null(prior@upper)

  collected <- matrix(0, nrow = 0, ncol = fit@estimator@n_dim_theta)
  n_tried <- 0L; batch <- 0L
  while (nrow(collected) < n && batch < max_sampling_batches) {
    batch <- batch + 1L
    draw_std <- de_sample(fit@estimator, xo_std, n)
    draw <- invert_standardizer(fit@std_theta, draw_std)
    n_tried <- n_tried + n
    if (bounded) draw <- draw[in_support(prior, draw), , drop = FALSE]
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
  out
}

method(log_prob, DirectPosterior) <- function(post, theta, x = NULL,
                                              normalize = TRUE,
                                              n_normalization = 10000L, ...) {
  fit <- post@fit
  theta <- as_theta_matrix(theta, fit@estimator@n_dim_theta)
  xo_std <- standardized_obs(post, x)
  theta_z <- apply_standardizer(fit@std_theta, theta)
  lp <- de_log_prob(fit@estimator, theta_z, xo_std) +
    standardizer_log_jac(fit@std_theta)

  prior <- fit@prior
  bounded <- !is.null(prior@lower) || !is.null(prior@upper)
  if (normalize && bounded) {
    draw_std <- de_sample(fit@estimator, xo_std, n_normalization)
    draw <- invert_standardizer(fit@std_theta, draw_std)
    acc <- max(mean(in_support(prior, draw)), 1 / n_normalization)
    lp <- lp - log(acc)
    lp[!in_support(prior, theta)] <- -Inf
  }
  lp
}

method(map_estimate, DirectPosterior) <- function(post, x = NULL,
                                                  n_init = 1000L, ...) {
  draws <- sample(post, n = n_init, obs = x)
  lp <- log_prob(post, draws, x = x, normalize = FALSE)
  start <- draws[which.max(lp), ]
  neg <- function(par) -log_prob(post, matrix(par, nrow = 1), x = x,
                                 normalize = FALSE)
  stats::optim(start, neg, method = "Nelder-Mead")$par
}

method(print, DirectPosterior) <- function(x, ...) {
  cat("<DirectPosterior>\n")
  cat(sprintf("  parameters (dim): %d\n", x@fit@estimator@n_dim_theta))
  cat(sprintf("  conditioned on x: %s\n",
              if (is.null(x@default_x)) "(none set)" else
                paste(signif(x@default_x[1, ], 4), collapse = ", ")))
  cat("  sample(post, n), log_prob(post, theta), map_estimate(post)\n")
  invisible(x)
}
