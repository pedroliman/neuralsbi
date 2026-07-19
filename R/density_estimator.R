#' Conditional density estimators
#'
#' A conditional density estimator learns \eqn{q_\phi(\theta \mid x)}. In
#' `neuralsbi` every estimator is trained in *standardized* space and exposes two
#' generics:
#'
#' * `de_log_prob(de, theta, x)` -- log density of `theta` given `x`
#' * `de_sample(de, x, n)` -- draw `n` parameter vectors given a single `x`
#'
#' Two estimators ship today:
#'
#' * `"mdn"` -- a Mixture Density Network (neural network -> Gaussian mixture),
#'   the workhorse, requires the `torch` back end.
#' * `"linear_gaussian"` -- a closed-form conditional Gaussian baseline
#'   (least-squares mean, residual covariance). No neural network, no `torch`.
#'   It is exact for linear-Gaussian simulators and doubles as a fast baseline
#'   and a regression-test oracle.
#'
#' @name density_estimator
NULL

#' @keywords internal
de_log_prob <- function(de, theta, x) UseMethod("de_log_prob")

#' @keywords internal
de_sample <- function(de, x, n) UseMethod("de_sample")

# ---- linear-Gaussian conditional estimator (pure R) -----------------------

#' @keywords internal
fit_linear_gaussian <- function(theta, x, ridge = 1e-6, verbose = FALSE) {
  theta <- as_theta_matrix(theta)
  x <- as_theta_matrix(x)
  n <- nrow(theta)
  p <- ncol(theta)
  X <- cbind(1, x)                       # design matrix with intercept
  # Ridge-regularized least squares: B = (X'X + rI)^-1 X'theta
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
  structure(
    list(B = B, Sigma = Sigma, chol = chol(Sigma), dim_theta = p),
    class = c("nsbi_de_lingauss", "nsbi_de")
  )
}

#' @keywords internal
lingauss_mean <- function(de, x) {
  x <- as_theta_matrix(x)
  cbind(1, x) %*% de$B
}

#' @export
de_log_prob.nsbi_de_lingauss <- function(de, theta, x) {
  theta <- as_theta_matrix(theta, de$dim_theta)
  mu <- lingauss_mean(de, x)
  if (nrow(mu) == 1L && nrow(theta) > 1L) {
    mu <- matrix(mu, nrow = nrow(theta), ncol = ncol(mu), byrow = TRUE)
  }
  dmvnorm_chol(theta, mu, de$chol, log = TRUE)
}

#' @export
de_sample.nsbi_de_lingauss <- function(de, x, n) {
  x <- as_theta_matrix(x)
  mu <- lingauss_mean(de, x)[1, ]
  z <- matrix(stats::rnorm(n * de$dim_theta), nrow = n)
  sweep(z %*% de$chol, 2, mu, `+`)
}

# ---- small multivariate-normal helpers ------------------------------------

#' Multivariate normal log density using a precomputed upper-Cholesky factor
#' (`R` such that `Sigma = t(R) %*% R`, i.e. `chol(Sigma)`).
#' @keywords internal
dmvnorm_chol <- function(x, mean, R, log = TRUE) {
  x <- as_theta_matrix(x)
  if (is.null(dim(mean))) mean <- matrix(mean, nrow = nrow(x), ncol = ncol(x),
                                         byrow = TRUE)
  d <- ncol(x)
  dev <- x - mean
  # Solve R' z = dev'  =>  quadratic form = colSums(z^2)
  z <- backsolve(R, t(dev), transpose = TRUE)
  quad <- colSums(z^2)
  logdet <- 2 * sum(log(diag(R)))
  out <- -0.5 * (d * log(2 * pi) + logdet + quad)
  if (log) out else exp(out)
}
