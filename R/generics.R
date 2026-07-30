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
#' @param ... Passed to the dispatched [sample()] method.
#' @return An `n x dim` matrix of posterior draws.
#' @export
sample_posterior <- function(post, n = 1000, obs = NULL, ...) {
  # Through the generic, not sample.nsbi_posterior() directly: an nle() fit's
  # posterior inherits nsbi_posterior but samples with MCMC, and its estimator
  # has the roles of theta and x swapped. The NPE sampler asked to draw from it
  # errors where dim_x and dim_theta differ, and where they happen to match it
  # returns draws from the wrong distribution without complaint.
  sample(post, n = n, obs = obs, ...)
}
