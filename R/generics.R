#' Draw samples
#'
#' `neuralsbi` promotes `sample()` to an S7 generic so that
#' `sample(posterior, n)` works in the style of the Python `sbi` package. Any
#' object without a dedicated method (vectors, etc.) falls back to
#' [base::sample()] unchanged.
#'
#' @param x Object to sample from (e.g. a [DirectPosterior]).
#' @param ... Passed to methods / [base::sample()].
#' @name sample
#' @return For a posterior, an `n x n_dim` matrix of draws; otherwise whatever
#'   [base::sample()] returns.
#' @export
method(sample, class_any) <- function(x, ...) base::sample(x, ...)

#' Sample from a posterior (non-generic alias)
#'
#' Identical to `sample(post, n)`; provided for users who prefer not to rely on
#' the generic.
#'
#' @param post A [DirectPosterior] object.
#' @param n Number of posterior draws.
#' @param obs Observation to condition on (defaults to the posterior's `x_obs`).
#' @param ... Passed to the [sample()] method.
#' @return An `n x n_dim` matrix of posterior draws.
#' @export
sample_posterior <- function(post, n = 1000, obs = NULL, ...) {
  sample(post, n = n, obs = obs, ...)
}
