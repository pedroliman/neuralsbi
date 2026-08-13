#' Summaries and tidy accessors
#'
#' `summary()` methods for fits, posteriors, and samples, plus
#' `as.data.frame()` for posterior draws so results drop straight into
#' data-frame workflows (dplyr, ggplot2, ...).
#'
#' @param object An `nsbi_samples` matrix from [sample()], an `nsbi_posterior`,
#'   or an `nsbi_npe`, `nsbi_nle` or `nsbi_nre` fit.
#' @param probs Quantiles to report.
#' @param n Number of draws used to summarize a posterior.
#' @param x Observation to condition on (defaults to the posterior's `x_obs`);
#'   for `as.data.frame()`, the samples object.
#' @param row.names,optional Standard [as.data.frame()] arguments.
#' @param ... Additional arguments passed to methods.
#' @return For samples and posteriors, a data frame with one row per parameter
#'   (mean, sd, and quantiles). For fits, an invisible list of training
#'   metadata. In that list `dim_x` counts the data dimension the estimator was
#'   trained on, which for an `nsbi_nle` or `nsbi_nre` fit is one observation
#'   rather than the whole data set.
#' @name summaries
NULL

#' @rdname summaries
#' @export
as.data.frame.nsbi_samples <- function(x, row.names = NULL, optional = FALSE, ...) {
  m <- unclass(x)
  attr(m, "acceptance_rate") <- NULL
  colnames(m) <- colnames(m) %||% paste0("theta", seq_len(ncol(m)))
  as.data.frame(m, row.names = row.names, optional = optional, ...)
}

#' @rdname summaries
#' @export
summary.nsbi_samples <- function(object,
                                 probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                                 ...) {
  m <- unclass(object)
  if (nrow(m) == 0L) {
    warning("summary(): 0 samples, so mean/sd/quantiles are all NA. ",
            "The posterior returned no draws inside the prior support; ",
            "see ?sample.nsbi_posterior.", call. = FALSE)
  }
  q <- t(apply(m, 2, stats::quantile, probs = probs))
  colnames(q) <- paste0("q", 100 * probs)
  out <- data.frame(
    parameter = colnames(m) %||% paste0("theta", seq_len(ncol(m))),
    mean = colMeans(m),
    sd = apply(m, 2, stats::sd),
    q, row.names = NULL, check.names = FALSE
  )
  out
}

#' @rdname summaries
#' @export
summary.nsbi_posterior <- function(object, n = 1000L, x = NULL, ...) {
  draws <- sample(object, n = n, obs = x)
  summary(draws, ...)
}

#' @rdname summaries
#' @export
summary.nsbi_npe <- function(object, ...) {
  fit_summary(object, list(density_estimator = object$density_estimator))
}

#' @rdname summaries
#' @export
# An NLE fit carries the same training metadata under the same names, and the
# print() call inside dispatches on the class, so there is nothing to
# specialize here.
summary.nsbi_nle <- function(object, ...) summary.nsbi_npe(object, ...)

#' @rdname summaries
#' @export
# An NRE fit records which classifier it trained rather than which density
# estimator, so that one field is named differently; everything else matches.
summary.nsbi_nre <- function(object, ...) {
  fit_summary(object, list(classifier = object$classifier))
}

#' Print a fit and return its training metadata invisibly
#'
#' The body every `summary()` method for a fit shares. `estimator` is the one
#' field whose name depends on what was trained: `density_estimator` for
#' [npe()] and [nle()], `classifier` for [nre()]. `print()` inside dispatches
#' on the object's own class, so nothing else here has to know which fit it has.
#'
#' @param object The fit.
#' @param estimator A one-element named list, spliced in ahead of the rest.
#' @keywords internal
fit_summary <- function(object, estimator) {
  info <- c(estimator, list(
    dim_theta = object$dim_theta,
    dim_x = object$dim_x,
    n_simulations = object$n_simulations,
    best_val_loss = object$de$best_val_loss %||% NA_real_,
    epochs_trained = if (!is.null(object$de$history)) nrow(object$de$history)
                     else NA_integer_
  ))
  print(object)
  if (!is.na(info$epochs_trained)) {
    cat(sprintf("  epochs trained    : %d\n", info$epochs_trained))
  }
  invisible(info)
}
