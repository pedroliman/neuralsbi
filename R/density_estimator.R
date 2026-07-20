#' Density estimators
#'
#' `neuralsbi` learns \eqn{q(\theta\mid x)} with an interchangeable estimator.
#' Two ship today, both trained in standardized space and implementing the
#' [de_log_prob()] / [de_sample()] generics:
#'
#' * [MDN] (`"mdn"`) -- neural Mixture Density Network (needs `torch`).
#' * [LinearGaussian] (`"linear_gaussian"`) -- closed-form conditional Gaussian;
#'   no neural network, no `torch`, exact for linear-Gaussian models. Doubles as
#'   a fast baseline and a regression-test oracle.
#'
#' @name density-estimators
NULL

# ---- linear-Gaussian conditional estimator (pure R) -----------------------

#' @keywords internal
fit_linear_gaussian <- function(theta, x, ridge = 1e-6, verbose = FALSE) {
  theta <- as_theta_matrix(theta)
  x <- as_theta_matrix(x)
  n <- nrow(theta)
  p <- ncol(theta)
  X <- cbind(1, x)                       # design matrix with intercept
  XtX <- crossprod(X)
  diag(XtX) <- diag(XtX) + ridge
  B <- solve(XtX, crossprod(X, theta))   # (q+1) x p
  mu <- X %*% B
  resid <- theta - mu
  Sigma <- crossprod(resid) / max(n - ncol(X), 1)
  diag(Sigma) <- diag(Sigma) + ridge
  verbose_cat(verbose, sprintf(
    "[linear_gaussian] fitted on %d sims, %d params, %d data dims\n",
    n, p, ncol(x)))
  LinearGaussian(n_dim_theta = p, n_dim_x = ncol(x),
                 B = B, Sigma = Sigma, chol = chol(Sigma))
}

#' @keywords internal
lingauss_mean <- function(de, x) {
  cbind(1, as_theta_matrix(x)) %*% de@B
}

method(de_log_prob, LinearGaussian) <- function(de, theta, x) {
  theta <- as_theta_matrix(theta, de@n_dim_theta)
  mu <- lingauss_mean(de, x)
  if (nrow(mu) == 1L && nrow(theta) > 1L) {
    mu <- matrix(mu, nrow = nrow(theta), ncol = ncol(mu), byrow = TRUE)
  }
  dmvnorm_chol(theta, mu, de@chol, log = TRUE)
}

method(de_sample, LinearGaussian) <- function(de, x, n) {
  mu <- lingauss_mean(de, x)[1, ]
  z <- matrix(stats::rnorm(n * de@n_dim_theta), nrow = n)
  sweep(z %*% de@chol, 2, mu, `+`)
}

# ---- small multivariate-normal helpers ------------------------------------

#' Multivariate normal log density using a precomputed upper-Cholesky factor
#' (`R` such that `Sigma = t(R) %*% R`, i.e. `chol(Sigma)`).
#' @keywords internal
dmvnorm_chol <- function(x, mean, R, log = TRUE) {
  x <- as_theta_matrix(x)
  if (is.null(dim(mean))) {
    mean <- matrix(mean, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
  }
  d <- ncol(x)
  dev <- x - mean
  z <- backsolve(R, t(dev), transpose = TRUE)
  quad <- colSums(z^2)
  logdet <- 2 * sum(log(diag(R)))
  out <- -0.5 * (d * log(2 * pi) + logdet + quad)
  if (log) out else exp(out)
}

method(print, DensityEstimator) <- function(x, ...) {
  cat(sprintf("<%s> theta_dim=%d, x_dim=%d\n",
              S7_class(x)@name, x@n_dim_theta, x@n_dim_x))
  invisible(x)
}
