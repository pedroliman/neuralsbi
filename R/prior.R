#' Priors for neural simulation-based inference
#'
#' A prior in `neuralsbi` is a lightweight object (class `nsbi_prior`) that knows
#' how to (a) draw samples and (b) evaluate its log-density. Bounded priors also
#' carry `lower`/`upper` support limits, which are used to reject out-of-support
#' posterior samples ("leakage" correction).
#'
#' @name priors
NULL

#' @keywords internal
new_prior <- function(sample_fn, log_prob_fn, dim, lower = NULL, upper = NULL,
                      type = "custom") {
  structure(
    list(
      sample = sample_fn,
      log_prob = log_prob_fn,
      dim = as.integer(dim),
      lower = lower,
      upper = upper,
      type = type
    ),
    class = "nsbi_prior"
  )
}

#' Box-uniform (independent uniform) prior
#'
#' @param low Numeric vector of lower bounds (one per parameter).
#' @param high Numeric vector of upper bounds (one per parameter).
#' @return An `nsbi_prior` object.
#' @examples
#' prior <- prior_uniform(low = c(-2, -2, -2), high = c(2, 2, 2))
#' theta <- sample_prior(prior, 5)
#' @export
prior_uniform <- function(low, high) {
  low <- as.numeric(low)
  high <- as.numeric(high)
  if (length(low) != length(high)) {
    stop("`low` and `high` must have the same length.", call. = FALSE)
  }
  if (any(high <= low)) {
    stop("Every `high` must be strictly greater than the matching `low`.",
         call. = FALSE)
  }
  d <- length(low)
  sample_fn <- function(n) {
    out <- matrix(stats::runif(n * d), nrow = n, ncol = d)
    sweep(sweep(out, 2, high - low, `*`), 2, low, `+`)
  }
  log_prob_fn <- function(theta) {
    theta <- as_theta_matrix(theta, d)
    inside <- rowSums(
      sweep(theta, 2, low, `>=`) & sweep(theta, 2, high, `<=`)
    ) == d
    const <- -sum(log(high - low))
    ifelse(inside, const, -Inf)
  }
  new_prior(sample_fn, log_prob_fn, d, lower = low, upper = high,
            type = "uniform")
}

#' Independent normal prior
#'
#' @param mean Numeric vector of means (one per parameter).
#' @param sd Numeric scalar or vector of standard deviations.
#' @return An `nsbi_prior` object.
#' @examples
#' prior <- prior_normal(mean = c(0, 0), sd = 1)
#' @export
prior_normal <- function(mean, sd = 1) {
  mean <- as.numeric(mean)
  d <- length(mean)
  sd <- as.numeric(sd)
  if (length(sd) == 1L) sd <- rep(sd, d)
  if (length(sd) != d) {
    stop("`sd` must be length 1 or the same length as `mean`.", call. = FALSE)
  }
  if (any(sd <= 0)) stop("`sd` must be positive.", call. = FALSE)
  sample_fn <- function(n) {
    z <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
    sweep(sweep(z, 2, sd, `*`), 2, mean, `+`)
  }
  log_prob_fn <- function(theta) {
    theta <- as_theta_matrix(theta, d)
    lp <- matrix(stats::dnorm(theta, mean = rep(mean, each = nrow(theta)),
                              sd = rep(sd, each = nrow(theta)), log = TRUE),
                 nrow = nrow(theta))
    rowSums(lp)
  }
  new_prior(sample_fn, log_prob_fn, d, lower = NULL, upper = NULL,
            type = "normal")
}

#' Build a prior from arbitrary sampling / density functions
#'
#' @param sample_fn Function `function(n)` returning an `n x dim` matrix.
#' @param log_prob_fn Function `function(theta)` returning a length-`n` vector of
#'   log densities. Optional; required only for methods/diagnostics that need it.
#' @param dim Number of parameters.
#' @param lower,upper Optional support bounds (numeric vectors) enabling
#'   out-of-support rejection.
#' @return An `nsbi_prior` object.
#' @export
prior_custom <- function(sample_fn, log_prob_fn = NULL, dim, lower = NULL,
                         upper = NULL) {
  if (is.null(log_prob_fn)) {
    log_prob_fn <- function(theta) rep(NA_real_, nrow(as_theta_matrix(theta, dim)))
  }
  new_prior(sample_fn, log_prob_fn, dim, lower = lower, upper = upper,
            type = "custom")
}

#' Draw samples from a prior
#'
#' @param prior An `nsbi_prior` object.
#' @param n Number of samples.
#' @return An `n x dim` matrix of parameter draws.
#' @export
sample_prior <- function(prior, n) {
  stopifnot(inherits(prior, "nsbi_prior"))
  out <- prior$sample(n)
  as_theta_matrix(out, prior$dim)
}

#' Test whether parameters lie within the prior support
#'
#' @param prior An `nsbi_prior` object.
#' @param theta A matrix (or vector) of parameters.
#' @return Logical vector, one entry per row of `theta`.
#' @export
within_support <- function(prior, theta) {
  theta <- as_theta_matrix(theta, prior$dim)
  if (is.null(prior$lower) && is.null(prior$upper)) {
    return(rep(TRUE, nrow(theta)))
  }
  ok <- rep(TRUE, nrow(theta))
  if (!is.null(prior$lower)) {
    ok <- ok & rowSums(sweep(theta, 2, prior$lower, `>=`)) == prior$dim
  }
  if (!is.null(prior$upper)) {
    ok <- ok & rowSums(sweep(theta, 2, prior$upper, `<=`)) == prior$dim
  }
  ok
}

#' @export
print.nsbi_prior <- function(x, ...) {
  cat(sprintf("<nsbi_prior> type=%s, dim=%d\n", x$type, x$dim))
  if (!is.null(x$lower)) {
    cat("  lower:", paste(signif(x$lower, 4), collapse = ", "), "\n")
    cat("  upper:", paste(signif(x$upper, 4), collapse = ", "), "\n")
  }
  invisible(x)
}
