#' Priors for neural simulation-based inference
#'
#' A prior in `neuralsbi` is a lightweight object (class `nsbi_prior`) that knows
#' how to (a) draw samples and (b) evaluate its log-density. Bounded priors also
#' carry `lower`/`upper` support limits, which are used to reject out-of-support
#' posterior samples ("leakage" correction).
#'
#' There are four ways to build one. [prior_uniform()] is the box prior most
#' benchmark tasks use. [prior_normal()] and the named families in
#' [prior_families] (log-normal, exponential, gamma, beta, Student-t, Cauchy
#' and the half versions of the last two) give one marginal per parameter,
#' under Stan's argument names. [prior_independent()] multiplies those together
#' into a joint prior, and [prior_truncated()] bounds one, renormalizing the
#' density by the mass it keeps. [prior_custom()] takes a sampler and a density
#' you write yourself, for anything the families do not cover.
#'
#' Prefer a named family over [prior_custom()] where one fits. A family carries
#' its support bounds into the leakage correction, survives [prior_truncated()]
#' with an exact normalizing constant, and is what [stan_code()] writes out as
#' a sampling statement; a custom prior does none of those.
#'
#' @seealso [prior_families], [prior_independent()], [prior_truncated()],
#'   [sample_prior()], [within_support()].
#' @name priors
NULL

#' @keywords internal
new_prior <- function(sample_fn, log_prob_fn, dim, lower = NULL, upper = NULL,
                      type = "custom", param_names = NULL, params = NULL) {
  structure(
    list(
      sample = sample_fn,
      log_prob = log_prob_fn,
      dim = as.integer(dim),
      lower = lower,
      upper = upper,
      type = type,
      param_names = param_names,
      # Distribution parameters kept alongside the closures, for code that has
      # to restate the prior somewhere else -- stan_code() writes it out as a
      # Stan sampling statement.
      params = params
    ),
    class = "nsbi_prior"
  )
}

#' Does this prior need bounding before it can be sampled or exported?
#'
#' [prior_uniform()] allows an infinite `low`/`high` bound so the result can
#' be handed to [prior_truncated()] later; until then it is an improper
#' distribution with no density over the whole real line. Two call sites need
#' to detect that and refuse it, each with a message suited to its own entry
#' point: the Stan export path (`check_finite_uniform_bounds()` in
#' `R/stan.R`) and the surrogate-potential path that backs MCMC sampling
#' (`surrogate_potential()` in `R/likelihood.R`). The check lives here once so
#' neither reimplements it.
#' @param prior An `nsbi_prior`.
#' @keywords internal
is_improper_uniform_prior <- function(prior) {
  identical(prior$type, "uniform") &&
    !(all(is.finite(prior$lower)) && all(is.finite(prior$upper)))
}

#' Box-uniform (independent uniform) prior
#'
#' @param low Numeric vector of lower bounds (one per parameter). Naming the
#'   vector (e.g. `c(beta = 0, gamma = 0)`) attaches those names to every
#'   downstream parameter matrix, posterior sample, and diagnostic plot.
#'   `NA`/`NaN` are rejected; `-Inf` is accepted (matching `high`'s `Inf`)
#'   but makes the prior improper, so [sample_prior()] on it errors unless
#'   it is first bounded with [prior_truncated()].
#' @param high Numeric vector of upper bounds (one per parameter). Same
#'   finiteness rule as `low`, with `Inf` in place of `-Inf`.
#' @return An `nsbi_prior` object.
#' @examples
#' prior <- prior_uniform(low = c(-2, -2, -2), high = c(2, 2, 2))
#' theta <- sample_prior(prior, 5)
#' @export
prior_uniform <- function(low, high) {
  param_names <- names(low) %||% names(high)
  check_finite(low, "low", allow_inf = TRUE)
  check_finite(high, "high", allow_inf = TRUE)
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
    if (any(!is.finite(high - low))) {
      stop(paste0("sample_prior() cannot draw from this prior directly: an ",
                  "infinite `low`/`high` bound makes it an improper ",
                  "distribution with no density to sample from. Wrap it in ",
                  "prior_truncated() to bound it to a finite interval first."),
           call. = FALSE)
    }
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
  # The closures above are the ones this prior has always used; `marginals` is
  # the same distribution restated in the canonical form of R/prior_families.R,
  # so prior_truncated() and prior_independent() can take a uniform apart.
  # An infinite bound is an improper prior with no CDF to invert, so it gets no
  # marginal; the closures above still work, as they always have.
  marginals <- if (all(is.finite(low)) && all(is.finite(high))) {
    lapply(seq_len(d), function(j) {
      new_marginal("uniform", list(min = low[[j]], max = high[[j]]),
                   label = param_label(param_names, j))
    })
  }
  new_prior(sample_fn, log_prob_fn, d, lower = low, upper = high,
            type = "uniform", param_names = param_names,
            params = list(low = low, high = high, marginals = marginals))
}

#' Independent normal prior
#'
#' @param mean Numeric vector of means (one per parameter). Naming the vector
#'   (e.g. `c(beta = 0, gamma = 0)`) attaches those names to every downstream
#'   parameter matrix, posterior sample, and diagnostic plot. Must be finite;
#'   a normal prior has no legitimate use for an infinite mean.
#' @param sd Numeric scalar or vector of standard deviations. Must be finite
#'   and positive.
#' @return An `nsbi_prior` object.
#' @examples
#' prior <- prior_normal(mean = c(0, 0), sd = 1)
#' @export
prior_normal <- function(mean, sd = 1) {
  param_names <- names(mean)
  check_finite(mean, "mean")
  mean <- as.numeric(mean)
  d <- length(mean)
  check_finite(sd, "sd")
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
  marginals <- if (all(is.finite(mean)) && all(is.finite(sd))) {
    lapply(seq_len(d), function(j) {
      new_marginal("normal", list(mean = mean[[j]], sd = sd[[j]]),
                   label = param_label(param_names, j))
    })
  }
  new_prior(sample_fn, log_prob_fn, d, lower = NULL, upper = NULL,
            type = "normal", param_names = param_names,
            params = list(mean = mean, sd = sd, marginals = marginals))
}

#' Build a prior from arbitrary sampling / density functions
#'
#' This is the one prior a user writes by hand, so it is the one where a
#' mistake is most likely, and the mistakes are quiet. A `log_prob_fn` that
#' returns a single number instead of one per row is legal R that surfaces much
#' later as an MCMC initialization failure; a `lower` of the wrong length is
#' recycled by `sweep()` into a support test that rejects the wrong draws.
#' Everything is therefore checked at construction, including one probe call of
#' `sample_fn(2)` and, when it is given, `log_prob_fn()` on those two rows.
#'
#' @param sample_fn Function `function(n)` returning an `n x dim` matrix, one
#'   row per draw.
#' @param log_prob_fn Function `function(theta)` returning a length-`n` vector of
#'   log densities, one per row of `theta`. Optional; required only for
#'   methods/diagnostics that need it.
#' @param dim Number of parameters.
#' @param lower,upper Optional support bounds enabling out-of-support
#'   rejection. Numeric of length `dim`, or length 1 to apply the same bound to
#'   every parameter.
#' @param param_names Optional character vector of parameter names, one per
#'   parameter. These name the columns of every downstream parameter matrix,
#'   posterior sample and diagnostic plot, and they decide how the simulator is
#'   called: a simulator whose formals are exactly these names receives one
#'   scalar per formal, and any other simulator receives the parameter vector as
#'   its first argument. See [nsbi_simulator].
#' @return An `nsbi_prior` object.
#' @section Stan export:
#' A custom prior is arbitrary R code rather than a named distribution with
#' parameters, so [stan_code()] cannot restate it as a Stan sampling statement.
#' Take `stan_code(fit, model = FALSE)` and write the model block yourself, or
#' build the prior out of [prior_families] and [prior_independent()] instead,
#' which [stan_code()] does write out.
#' @examples
#' prior <- prior_custom(
#'   sample_fn = function(n) cbind(rexp(n, 1), rexp(n, 2)),
#'   log_prob_fn = function(theta) {
#'     dexp(theta[, 1], 1, log = TRUE) + dexp(theta[, 2], 2, log = TRUE)
#'   },
#'   dim = 2, lower = 0, param_names = c("beta", "gamma")
#' )
#' @export
prior_custom <- function(sample_fn, log_prob_fn = NULL, dim, lower = NULL,
                         upper = NULL, param_names = NULL) {
  d <- check_count(dim, "dim", min = 1L, why = "(the number of parameters)")
  check_function(sample_fn, "sample_fn", what = "the number of draws")
  if (!is.null(log_prob_fn)) {
    check_function(log_prob_fn, "log_prob_fn", what = "a matrix of parameters")
  }
  lower <- check_bound(lower, "lower", d)
  upper <- check_bound(upper, "upper", d)
  if (!is.null(lower) && !is.null(upper) && any(upper <= lower)) {
    stop("Every `upper` must be strictly greater than the matching `lower`.",
         call. = FALSE)
  }
  ok_names <- is.null(param_names) ||
    (is.character(param_names) && length(param_names) == d &&
       !anyNA(param_names) && all(nzchar(param_names)))
  if (!ok_names) {
    stop(sprintf(paste0("`param_names` must be a character vector with one ",
                        "non-empty name per parameter (%s), not %s."),
                 n_things(d, "parameter"), describe_value(param_names)),
         call. = FALSE)
  }

  # One probe call, at construction rather than per draw. Whatever comes back
  # has to survive as_theta_matrix() in sample_prior(), so coerce it the same
  # way and check the shape here, where the error can still name `sample_fn`.
  probe <- tryCatch(sample_fn(2L), error = function(e) {
    stop(sprintf("`sample_fn` failed when called as `sample_fn(2)`: %s",
                 conditionMessage(e)),
         call. = FALSE)
  })
  probe <- as_theta_matrix(check_numeric(probe, "sample_fn(2)"))
  if (!identical(dim(probe), c(2L, d))) {
    stop(sprintf(paste0("`sample_fn(2)` must return a 2 x %d matrix (one row ",
                        "per draw, one parameter per column), but it returned ",
                        "a %d x %d matrix. See ?prior_custom."),
                 d, nrow(probe), ncol(probe)),
         call. = FALSE)
  }

  if (is.null(log_prob_fn)) {
    log_prob_fn <- function(theta) rep(NA_real_, nrow(as_theta_matrix(theta, d)))
  } else {
    lp <- tryCatch(log_prob_fn(probe), error = function(e) {
      stop(sprintf(paste0("`log_prob_fn` failed on a 2 x %d matrix of draws ",
                          "from `sample_fn`: %s"), d, conditionMessage(e)),
           call. = FALSE)
    })
    if (!is.numeric(lp) || length(lp) != 2L) {
      got <- if (is.null(lp)) "NULL" else
        sprintf("a length-%d %s vector", length(lp), class(lp)[1L])
      stop(sprintf(paste0("`log_prob_fn` must return one log-density per row ",
                          "of `theta`, but on a 2-row matrix it returned %s. ",
                          "See ?prior_custom."), got),
           call. = FALSE)
    }
  }

  new_prior(sample_fn, log_prob_fn, d, lower = lower, upper = upper,
            type = "custom", param_names = param_names)
}

#' Draw samples from a prior
#'
#' @param prior An `nsbi_prior` object.
#' @param n Number of samples.
#' @return An `n x dim` matrix of parameter draws.
#' @export
sample_prior <- function(prior, n) {
  stopifnot(inherits(prior, "nsbi_prior"))
  n <- check_count(n, "n")
  out <- as_theta_matrix(prior$sample(n), prior$dim)
  if (is.null(colnames(out)) && !is.null(prior$param_names)) {
    colnames(out) <- prior$param_names
  }
  out
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
  if (!is.null(x$param_names)) {
    cat("  parameters:", paste(x$param_names, collapse = ", "), "\n")
  }
  # Printed separately: a custom prior can bound one side and not the other,
  # and within_support() has always allowed that.
  if (!is.null(x$lower)) {
    cat("  lower:", paste(signif(x$lower, 4), collapse = ", "), "\n")
  }
  if (!is.null(x$upper)) {
    cat("  upper:", paste(signif(x$upper, 4), collapse = ", "), "\n")
  }
  # One line per marginal, which is the only place a family's parameters show
  # up. Skipped for uniform and normal: a uniform is fully described by the
  # bounds already printed, and prior_normal() has printed the way it does
  # since before families existed.
  if (!is.null(x$params$marginals) && !x$type %in% c("uniform", "normal")) {
    described <- describe_marginals(x$params$marginals)
    nm <- x$param_names %||% sprintf("theta[%d]", seq_along(described))
    for (j in seq_along(described)) {
      cat(sprintf("  %s ~ %s\n", nm[[j]], described[[j]]))
    }
  }
  invisible(x)
}
