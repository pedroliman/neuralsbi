#' Standardization (z-scoring) helpers
#'
#' Neural density estimators train far more reliably when inputs and targets are
#' standardized to roughly zero mean and unit variance. `neuralsbi` learns these
#' transforms from the training simulations, applies them internally, and inverts
#' them when returning posterior draws / densities.
#'
#' @name standardize
#' @keywords internal
NULL

#' Learn a standardizer from a matrix
#'
#' A column with no spread cannot be divided by its standard deviation, so it
#' keeps scale 1. Pass `what` to hear about it: the guard is silent otherwise,
#' and a constant column is worth a word because nothing downstream will
#' complain. Training converges, the posterior looks plausible, and the
#' coordinate does nothing.
#'
#' @param x The matrix to learn from.
#' @param eps Standard deviations below this count as no spread.
#' @param what Name of the argument `x` came from (`"theta"` or `"x"`), used in
#'   the warning. `NULL`, the default, warns about nothing. The
#'   `standardize = FALSE` path in [prepare_simulations()] builds a degenerate
#'   standardizer from a one-row zero matrix on purpose, and the diagnostics
#'   standardize draws they generated themselves.
#' @keywords internal
fit_standardizer <- function(x, eps = 1e-8, what = NULL) {
  x <- as_theta_matrix(x)
  center <- colMeans(x)
  scale <- apply(x, 2, stats::sd)
  flat <- scale < eps | !is.finite(scale)
  scale[flat] <- 1
  if (!is.null(what) && any(flat)) warn_constant_columns(x, flat, what)
  structure(list(center = center, scale = scale), class = "nsbi_standardizer")
}

#' Warn about columns standardization cannot scale
#'
#' Names the columns by name where they have one and by index otherwise, since
#' an index is no help when the matrix came from a data frame and a name is not
#' available when it did not. This is a warning rather than an error because a
#' single row is a legitimate way to get here: `sd()` of one value is `NA`.
#'
#' @param x The matrix being standardized.
#' @param flat Logical vector marking the columns with no usable spread.
#' @param what Name of the argument `x` came from, `"theta"` or `"x"`.
#' @keywords internal
warn_constant_columns <- function(x, flat, what) {
  nm <- colnames(x)
  labels <- vapply(which(flat), function(j) {
    if (!is.null(nm) && !is.na(nm[j]) && nzchar(nm[j])) nm[j] else
      sprintf("column %d", j)
  }, character(1))
  shown <- paste(labels, collapse = ", ")
  msg <- if (nrow(x) < 2L) {
    sprintf(paste0("`%s` has one row, so the standard deviation of %s (%s) is ",
                   "undefined. They are left unscaled; see ?standardize."),
            what, n_things(sum(flat), "column"), shown)
  } else {
    advice <- if (identical(what, "x")) {
      paste("A summary statistic that does not vary carries no information",
            "about `theta`, and the estimator will spend capacity on it.")
    } else {
      paste("A parameter that does not vary across the simulations cannot be",
            "identified from them.")
    }
    sprintf("`%s` has %s (%s) holding the same value in every row. %s See ?standardize.",
            what, n_things(sum(flat), "constant column"), shown, advice)
  }
  warning(msg, call. = FALSE)
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
