#' Standardization (z-scoring) helpers
#'
#' Neural density estimators train far more reliably when inputs and targets are
#' standardized to roughly zero mean and unit variance. `neuralsbi` learns these
#' transforms from the training simulations, applies them internally, and inverts
#' them when returning posterior draws / densities. This mirrors what the Python
#' `sbi` package does with its z-scoring transforms.
#'
#' @name standardize
#' @keywords internal
NULL

#' @keywords internal
fit_standardizer <- function(x, eps = 1e-8) {
  x <- as_theta_matrix(x)
  center <- colMeans(x)
  scale <- apply(x, 2, stats::sd)
  scale[scale < eps | !is.finite(scale)] <- 1
  structure(list(center = center, scale = scale), class = "nsbi_standardizer")
}

#' @keywords internal
apply_standardizer <- function(std, x) {
  x <- as_theta_matrix(x, length(std$center))
  sweep(sweep(x, 2, std$center, `-`), 2, std$scale, `/`)
}

#' @keywords internal
invert_standardizer <- function(std, z) {
  z <- as_theta_matrix(z, length(std$center))
  sweep(sweep(z, 2, std$scale, `*`), 2, std$center, `+`)
}

#' Log absolute Jacobian determinant of the *inverse* standardization
#' (standardized -> original). Constant, so a scalar.
#' @keywords internal
standardizer_log_jac <- function(std) {
  -sum(log(std$scale))
}
