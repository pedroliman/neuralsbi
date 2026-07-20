#' Box-uniform (independent uniform) prior
#'
#' @param low Numeric vector of lower bounds (one per parameter).
#' @param high Numeric vector of upper bounds (one per parameter).
#' @return A [UniformPrior] object.
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
  UniformPrior(n_dim = length(low), lower = low, upper = high,
               low = low, high = high)
}

#' Independent normal prior
#'
#' @param mean Numeric vector of means (one per parameter).
#' @param sd Numeric scalar or vector of standard deviations.
#' @return A [NormalPrior] object.
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
  NormalPrior(n_dim = d, mean = mean, sd = sd)
}

#' Build a prior from arbitrary sampling / density functions
#'
#' @param sample_fn Function `function(n)` returning an `n x dim` matrix.
#' @param log_prob_fn Function `function(theta)` returning a length-`n` vector of
#'   log densities. Optional; required only for methods/diagnostics that need it.
#' @param dim Number of parameters.
#' @param lower,upper Optional support bounds (numeric vectors) enabling
#'   out-of-support rejection.
#' @return A [CustomPrior] object.
#' @export
prior_custom <- function(sample_fn, log_prob_fn = NULL, dim, lower = NULL,
                         upper = NULL) {
  if (is.null(log_prob_fn)) {
    log_prob_fn <- function(theta) {
      rep(NA_real_, nrow(as_theta_matrix(theta, dim)))
    }
  }
  CustomPrior(n_dim = as.integer(dim), lower = lower, upper = upper,
              sample_fn = sample_fn, log_prob_fn = log_prob_fn)
}

# ---- draw_prior -----------------------------------------------------------

method(draw_prior, UniformPrior) <- function(prior, n) {
  d <- prior@n_dim
  out <- matrix(stats::runif(n * d), nrow = n, ncol = d)
  sweep(sweep(out, 2, prior@high - prior@low, `*`), 2, prior@low, `+`)
}

method(draw_prior, NormalPrior) <- function(prior, n) {
  d <- prior@n_dim
  z <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
  sweep(sweep(z, 2, prior@sd, `*`), 2, prior@mean, `+`)
}

method(draw_prior, CustomPrior) <- function(prior, n) {
  as_theta_matrix(prior@sample_fn(n), prior@n_dim)
}

# ---- prior_log_prob -------------------------------------------------------

method(prior_log_prob, UniformPrior) <- function(prior, theta) {
  theta <- as_theta_matrix(theta, prior@n_dim)
  inside <- rowSums(sweep(theta, 2, prior@low, `>=`) &
                      sweep(theta, 2, prior@high, `<=`)) == prior@n_dim
  const <- -sum(log(prior@high - prior@low))
  ifelse(inside, const, -Inf)
}

method(prior_log_prob, NormalPrior) <- function(prior, theta) {
  theta <- as_theta_matrix(theta, prior@n_dim)
  lp <- matrix(stats::dnorm(theta,
                            mean = rep(prior@mean, each = nrow(theta)),
                            sd = rep(prior@sd, each = nrow(theta)),
                            log = TRUE), nrow = nrow(theta))
  rowSums(lp)
}

method(prior_log_prob, CustomPrior) <- function(prior, theta) {
  prior@log_prob_fn(as_theta_matrix(theta, prior@n_dim))
}

# ---- in_support (shared, uses box bounds) ---------------------------------

method(in_support, Prior) <- function(prior, theta) {
  theta <- as_theta_matrix(theta, prior@n_dim)
  if (is.null(prior@lower) && is.null(prior@upper)) {
    return(rep(TRUE, nrow(theta)))
  }
  ok <- rep(TRUE, nrow(theta))
  if (!is.null(prior@lower)) {
    ok <- ok & rowSums(sweep(theta, 2, prior@lower, `>=`)) == prior@n_dim
  }
  if (!is.null(prior@upper)) {
    ok <- ok & rowSums(sweep(theta, 2, prior@upper, `<=`)) == prior@n_dim
  }
  ok
}

# ---- user-facing convenience wrappers -------------------------------------

#' Draw samples from a prior
#'
#' @param prior A [Prior] object.
#' @param n Number of samples.
#' @return An `n x n_dim` matrix of parameter draws.
#' @export
sample_prior <- function(prior, n) {
  as_theta_matrix(draw_prior(prior, n), prior@n_dim)
}

#' Test whether parameters lie within the prior support
#'
#' @param prior A [Prior] object.
#' @param theta A matrix (or vector) of parameters.
#' @return Logical vector, one entry per row of `theta`.
#' @export
within_support <- function(prior, theta) in_support(prior, theta)

method(print, Prior) <- function(x, ...) {
  cat(sprintf("<%s> n_dim=%d\n", S7_class(x)@name, x@n_dim))
  if (!is.null(x@lower)) {
    cat("  lower:", paste(signif(x@lower, 4), collapse = ", "), "\n")
    cat("  upper:", paste(signif(x@upper, 4), collapse = ", "), "\n")
  }
  invisible(x)
}
