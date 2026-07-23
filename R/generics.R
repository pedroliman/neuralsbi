#' Draw samples (S3 generic)
#'
#' `neuralsbi` turns [base::sample()] into an S3 generic so that
#' `sample(posterior, n)` reads the way statisticians expect. For
#' any object without a dedicated method (vectors, etc.) this falls back to
#' [base::sample()] unchanged.
#'
#' @param x Object to sample from.
#' @param ... Passed on to methods / [base::sample()].
#' @return Whatever the dispatched method returns. The default method returns
#'   the result of [base::sample()]; [sample.nsbi_posterior()] returns an
#'   `n x dim` matrix of posterior draws.
#' @export
sample <- function(x, ...) UseMethod("sample")

#' @rdname sample
#' @export
sample.default <- function(x, ...) base::sample(x, ...)

#' Sample from a posterior (non-generic alias)
#'
#' Identical to `sample(post, n)`; provided for users who prefer not to rely on
#' the generic.
#'
#' @param post An `nsbi_posterior` object.
#' @param n Number of posterior draws.
#' @param obs Observation to condition on (defaults to the posterior's `x_obs`).
#' @param ... Passed to [sample.nsbi_posterior()].
#' @return An `n x dim` matrix of posterior draws.
#' @export
sample_posterior <- function(post, n = 1000, obs = NULL, ...) {
  sample.nsbi_posterior(post, n = n, obs = obs, ...)
}
