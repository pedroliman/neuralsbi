#' Named prior families
#'
#' Priors built from a named distribution family, one marginal per parameter,
#' following Stan's argument names and argument order. Each constructor is
#' vectorized: pass a vector where Stan would write one sampling statement per
#' component, and pass a scalar to use the same value for every parameter.
#' Naming any of the vectors (e.g. `c(beta = 2, gamma = 3)`) names the
#' parameters, exactly as naming `low` does in [prior_uniform()].
#'
#' Families with constrained support set `lower`/`upper` on the prior, so the
#' posterior's leakage correction (see [posterior()]) rejects and renormalizes
#' against the right region without any further declaration. A gamma prior
#' bounds at zero, a beta prior at zero and one, and the half families at zero.
#'
#' A prior built here can be truncated further with [prior_truncated()],
#' combined with others through [prior_independent()], and written out by
#' [stan_code()] as the sampling statement it came from. That is the reason
#' these exist as families rather than as recipes for [prior_custom()]: the
#' family is what carries over to Stan and to the truncation constant.
#'
#' `prior_half_normal()` and `prior_half_cauchy()` are the zero-centred
#' families truncated to \eqn{[0, \infty)}, matching Stan's idiom of declaring
#' `real<lower=0>` and writing `sigma ~ normal(0, 1)`. Their densities carry
#' the `log(2)` renormalization, which Stan drops as a constant and this
#' package cannot. For a half-normal centred somewhere other than zero, wrap
#' [prior_normal()] in [prior_truncated()].
#'
#' @param meanlog,sdlog Log-scale mean and standard deviation of the
#'   log-normal, as in [stats::dlnorm()] and Stan's `lognormal`.
#' @param rate Rate of the exponential or gamma, as in [stats::dexp()] and
#'   Stan's `exponential`/`gamma`. Must be positive.
#' @param shape Shape of the gamma. Must be positive.
#' @param shape1,shape2 Beta shape parameters, as in [stats::dbeta()]. Both
#'   must be positive.
#' @param df Degrees of freedom of the Student-t. Must be positive.
#' @param location,scale Location and scale of the Student-t, the Cauchy or the
#'   half-Cauchy. `scale` must be positive.
#' @param sd Standard deviation of the normal underlying the half-normal. Must
#'   be positive.
#' @return An `nsbi_prior` object.
#' @seealso [prior_independent()] to combine several of these into one joint
#'   prior, [prior_truncated()] to bound one, and [priors] for the rest.
#' @examples
#' # One log-normal per parameter, named through the vector.
#' prior <- prior_lognormal(meanlog = c(beta = log(0.4), gamma = log(0.125)),
#'                          sdlog = c(0.5, 0.2))
#' theta <- sample_prior(prior, 5)
#'
#' # A scalar argument is shared: three gammas with the same rate.
#' prior_gamma(shape = c(2, 5, 9), rate = 3)
#'
#' # Support bounds come from the family, and drive the leakage correction.
#' within_support(prior_beta(2, 15), c(0.1, 1.5))
#' @name prior_families
NULL

# ---- the family registry --------------------------------------------------

#' The distribution families a prior can be built from
#'
#' One entry per family: the parameter names in Stan's order, the R
#' density/CDF/quantile trio behind them, and the natural support. Everything
#' else in this file is generic over that entry, which is why registering a
#' family is all it takes for [prior_truncated()] to renormalize it and for
#' [stan_code()] to write it out.
#'
#' The `d`, `p` and `q` functions are always called with the value positionally
#' and the distribution parameters by name, since R spells the first argument
#' `x`, `q` and `p` in turn while the parameters keep their names throughout.
#'
#' @param family Family name, or anything else to get `NULL` back.
#' @return A list, or `NULL` for an unregistered name.
#' @keywords internal
prior_family <- function(family) {
  switch(
    family,
    uniform = list(
      args = c("min", "max"), stan = "uniform",
      d = stats::dunif, p = stats::punif, q = stats::qunif,
      support = function(a) c(a$min, a$max)
    ),
    normal = list(
      args = c("mean", "sd"), stan = "normal",
      d = stats::dnorm, p = stats::pnorm, q = stats::qnorm,
      support = function(a) c(-Inf, Inf)
    ),
    lognormal = list(
      args = c("meanlog", "sdlog"), stan = "lognormal",
      d = stats::dlnorm, p = stats::plnorm, q = stats::qlnorm,
      support = function(a) c(0, Inf)
    ),
    exponential = list(
      args = "rate", stan = "exponential",
      d = stats::dexp, p = stats::pexp, q = stats::qexp,
      support = function(a) c(0, Inf)
    ),
    gamma = list(
      args = c("shape", "rate"), stan = "gamma",
      d = stats::dgamma, p = stats::pgamma, q = stats::qgamma,
      support = function(a) c(0, Inf)
    ),
    beta = list(
      args = c("shape1", "shape2"), stan = "beta",
      d = stats::dbeta, p = stats::pbeta, q = stats::qbeta,
      support = function(a) c(0, 1)
    ),
    cauchy = list(
      args = c("location", "scale"), stan = "cauchy",
      d = stats::dcauchy, p = stats::pcauchy, q = stats::qcauchy,
      support = function(a) c(-Inf, Inf)
    ),
    student_t = list(
      args = c("df", "location", "scale"), stan = "student_t",
      d = dstudent_t, p = pstudent_t, q = qstudent_t,
      support = function(a) c(-Inf, Inf)
    ),
    NULL
  )
}

#' Location-scale Student-t
#'
#' `stats::dt()` only knows the standardized t, and Stan's `student_t` is the
#' location-scale one. These three wrap the shift and the scale, including the
#' `-log(scale)` Jacobian on the density, and take their arguments under the
#' names the registry uses.
#' @param x,q,p Value, quantile or probability.
#' @param df,location,scale Distribution parameters.
#' @param log Return the log density.
#' @keywords internal
dstudent_t <- function(x, df, location = 0, scale = 1, log = FALSE) {
  out <- stats::dt((x - location) / scale, df = df, log = TRUE) - log(scale)
  if (isTRUE(log)) out else exp(out)
}

#' @rdname dstudent_t
#' @keywords internal
pstudent_t <- function(q, df, location = 0, scale = 1) {
  stats::pt((q - location) / scale, df = df)
}

#' @rdname dstudent_t
#' @keywords internal
qstudent_t <- function(p, df, location = 0, scale = 1) {
  location + scale * stats::qt(p, df = df)
}

# ---- one marginal ---------------------------------------------------------

#' One parameter's marginal: a family, its parameters, and its bounds
#'
#' The canonical form every named-family prior is stored in, under
#' `prior$params$marginals`. `lower`/`upper` are the effective bounds, that is
#' the family's own support intersected with whatever truncation was asked for,
#' and `log_norm` is the log of the probability mass left inside them. Holding
#' the truncated mass here rather than recomputing it per call is what keeps
#' the density proper without paying for a CDF evaluation on every row.
#'
#' @param family Family name, a key of [prior_family()].
#' @param args Named list of scalar distribution parameters.
#' @param lower,upper Truncation bounds, before intersecting with the family's
#'   own support.
#' @param label How to refer to this parameter in an error message.
#' @return A marginal specification.
#' @keywords internal
new_marginal <- function(family, args, lower = -Inf, upper = Inf,
                         label = "the parameter") {
  f <- prior_family(family)
  args <- args[f$args]
  natural <- f$support(args)
  lower <- max(lower, natural[1L])
  upper <- min(upper, natural[2L])
  if (upper <= lower) {
    stop(sprintf(paste0("The bounds for %s leave nothing of the %s ",
                        "distribution: [%s, %s] is empty once it is ",
                        "intersected with the family's support [%s, %s]."),
                 label, family, format(lower), format(upper),
                 format(natural[1L]), format(natural[2L])),
         call. = FALSE)
  }
  p_lower <- do.call(f$p, c(list(lower), args))
  p_upper <- do.call(f$p, c(list(upper), args))
  mass <- p_upper - p_lower
  if (!is.finite(mass) || mass <= 0) {
    stop(sprintf(paste0("The bounds for %s leave no probability mass: a %s ",
                        "prior puts %s of its mass on [%s, %s], so the ",
                        "truncated density has no normalizing constant."),
                 label, family, format(mass), format(lower), format(upper)),
         call. = FALSE)
  }
  list(family = family, args = args, lower = lower, upper = upper,
       natural = natural, p_lower = p_lower, p_upper = p_upper,
       log_norm = log(mass))
}

#' Log density of one marginal, renormalized and masked outside its bounds
#' @param m A marginal from [new_marginal()].
#' @param x Numeric vector of values.
#' @keywords internal
marginal_log_prob <- function(m, x) {
  f <- prior_family(m$family)
  out <- do.call(f$d, c(list(x), m$args, list(log = TRUE))) - m$log_norm
  # which() rather than a logical index so an NA value stays NA instead of
  # erroring on a missing subscript.
  out[which(x < m$lower | x > m$upper)] <- -Inf
  out
}

#' Draw from one marginal by inverting its CDF between the bounds
#'
#' Inverse-CDF sampling covers the truncated and untruncated cases with the
#' same line, and it is exact: no rejection loop that stalls when the bounds
#' cut into a tail.
#' @param m A marginal from [new_marginal()].
#' @param n Number of draws.
#' @keywords internal
marginal_sample <- function(m, n) {
  f <- prior_family(m$family)
  u <- m$p_lower + stats::runif(n) * (m$p_upper - m$p_lower)
  do.call(f$q, c(list(u), m$args))
}

#' Assemble marginals into an independent joint prior
#'
#' @param marginals List of marginals, one per parameter.
#' @param param_names Optional parameter names.
#' @param type The `type` recorded on the prior.
#' @return An `nsbi_prior` object.
#' @keywords internal
marginal_prior <- function(marginals, param_names = NULL,
                           type = "independent") {
  d <- length(marginals)
  marginals <- unname(marginals)
  lower <- unname(vapply(marginals, `[[`, numeric(1), "lower"))
  upper <- unname(vapply(marginals, `[[`, numeric(1), "upper"))

  sample_fn <- function(n) {
    out <- matrix(0, nrow = n, ncol = d)
    for (j in seq_len(d)) out[, j] <- marginal_sample(marginals[[j]], n)
    out
  }
  log_prob_fn <- function(theta) {
    theta <- as_theta_matrix(theta, d)
    lp <- matrix(0, nrow = nrow(theta), ncol = d)
    for (j in seq_len(d)) lp[, j] <- marginal_log_prob(marginals[[j]], theta[, j])
    rowSums(lp)
  }

  new_prior(sample_fn, log_prob_fn, d,
            lower = if (all(lower == -Inf)) NULL else lower,
            upper = if (all(upper == Inf)) NULL else upper,
            type = type, param_names = param_names,
            params = list(marginals = marginals))
}

# ---- argument handling shared by the constructors -------------------------

#' Recycle a family's arguments to one value per parameter
#'
#' `prior_gamma(shape = c(2, 3), rate = 1)` is two parameters with a shared
#' rate, and `prior_gamma(shape = c(2, 3), rate = c(1, 2, 3))` is a mistake.
#' The longest argument sets the number of parameters and every other one has
#' to be that length or length 1, which is `prior_normal()`'s rule for `sd`
#' generalized to families with more than one parameter.
#'
#' @param args Named list of numeric vectors.
#' @param fname Constructor name, for the error message.
#' @return `args`, each entry a length-`d` double vector without names.
#' @keywords internal
recycle_family_args <- function(args, fname) {
  lens <- lengths(args)
  d <- max(lens)
  bad <- names(args)[lens != 1L & lens != d]
  if (length(bad) > 0L) {
    stop(sprintf(paste0("Every argument to %s() must be length 1 or the ",
                        "length of the longest one (%d), but %s %s not."),
                 fname, d, paste(sprintf("`%s`", bad), collapse = " and "),
                 if (length(bad) == 1L) "is" else "are"),
         call. = FALSE)
  }
  lapply(args, function(a) rep_len(as.double(unname(a)), d))
}

#' Parameter names taken from whichever argument carries them
#'
#' `prior_uniform()` reads them off `low` or `high`; a family with three
#' arguments has three places to look, and the first complete set wins.
#' @param args Named list of the raw (un-recycled) arguments.
#' @param d Number of parameters.
#' @keywords internal
family_param_names <- function(args, d) {
  for (a in args) {
    nm <- names(a)
    if (!is.null(nm) && length(nm) == d && !anyNA(nm) && all(nzchar(nm))) {
      return(nm)
    }
  }
  NULL
}

#' Build a prior from one named family
#'
#' The body every constructor in this file shares: validate, recycle, turn each
#' parameter into a marginal, hand the lot to [marginal_prior()].
#'
#' @param family Family name, a key of [prior_family()].
#' @param args Named list of the validated arguments, in any order.
#' @param fname Constructor name, for error messages.
#' @param lower,upper Truncation applied to every parameter, for the half
#'   families.
#' @param type The `type` recorded on the prior; defaults to the family name.
#' @keywords internal
family_prior <- function(family, args, fname, lower = -Inf, upper = Inf,
                         type = family) {
  d <- max(lengths(args))
  param_names <- family_param_names(args, d)
  args <- recycle_family_args(args, fname)
  marginals <- lapply(seq_len(d), function(j) {
    new_marginal(family, lapply(args, `[[`, j), lower = lower, upper = upper,
                 label = param_label(param_names, j))
  })
  marginal_prior(marginals, param_names = param_names, type = type)
}

#' Name one parameter in an error message, by name where there is one
#' @param param_names Parameter names, or `NULL`.
#' @param j Parameter index.
#' @keywords internal
param_label <- function(param_names, j) {
  if (is.null(param_names)) sprintf("parameter %d", j) else
    sprintf("`%s`", param_names[[j]])
}

# ---- the constructors -----------------------------------------------------

#' @rdname prior_families
#' @export
prior_lognormal <- function(meanlog = 0, sdlog = 1) {
  family_prior(
    "lognormal",
    list(meanlog = check_family_param(meanlog, "meanlog"),
         sdlog = check_family_param(sdlog, "sdlog", positive = TRUE)),
    fname = "prior_lognormal"
  )
}

#' @rdname prior_families
#' @export
prior_exponential <- function(rate = 1) {
  family_prior(
    "exponential",
    list(rate = check_family_param(rate, "rate", positive = TRUE)),
    fname = "prior_exponential"
  )
}

#' @rdname prior_families
#' @export
prior_gamma <- function(shape, rate = 1) {
  family_prior(
    "gamma",
    list(shape = check_family_param(shape, "shape", positive = TRUE),
         rate = check_family_param(rate, "rate", positive = TRUE)),
    fname = "prior_gamma"
  )
}

#' @rdname prior_families
#' @export
prior_beta <- function(shape1, shape2) {
  family_prior(
    "beta",
    list(shape1 = check_family_param(shape1, "shape1", positive = TRUE),
         shape2 = check_family_param(shape2, "shape2", positive = TRUE)),
    fname = "prior_beta"
  )
}

#' @rdname prior_families
#' @export
prior_student_t <- function(df, location = 0, scale = 1) {
  family_prior(
    "student_t",
    list(df = check_family_param(df, "df", positive = TRUE),
         location = check_family_param(location, "location"),
         scale = check_family_param(scale, "scale", positive = TRUE)),
    fname = "prior_student_t"
  )
}

#' @rdname prior_families
#' @export
prior_cauchy <- function(location = 0, scale = 1) {
  family_prior(
    "cauchy",
    list(location = check_family_param(location, "location"),
         scale = check_family_param(scale, "scale", positive = TRUE)),
    fname = "prior_cauchy"
  )
}

#' @rdname prior_families
#' @export
prior_half_normal <- function(sd = 1) {
  family_prior(
    "normal",
    list(mean = 0, sd = check_family_param(sd, "sd", positive = TRUE)),
    fname = "prior_half_normal", lower = 0, type = "half_normal"
  )
}

#' @rdname prior_families
#' @export
prior_half_cauchy <- function(scale = 1) {
  family_prior(
    "cauchy",
    list(location = 0,
         scale = check_family_param(scale, "scale", positive = TRUE)),
    fname = "prior_half_cauchy", lower = 0, type = "half_cauchy"
  )
}

# ---- composition and truncation -------------------------------------------

#' Combine independent priors into one joint prior
#'
#' Stan writes a joint prior as one sampling statement per parameter and lets
#' the product take care of itself. This is that, as an object: give it the
#' per-parameter priors in the order the simulator expects them and it returns
#' the product prior, with `dim` the total, `lower`/`upper` stacked, and a
#' `log_prob` that sums the parts. It is the thing most [prior_custom()] calls
#' were written to do by hand.
#'
#' Components may themselves cover several parameters, so
#' `prior_independent(prior_normal(mean = c(0, 0)), prior_beta(2, 15))` is a
#' three-parameter prior. Any `nsbi_prior` works as a component, including a
#' [prior_custom()]; the result can only be written out by [stan_code()] when
#' every component comes from a named family (see [prior_families]).
#'
#' @param ... Priors, one per block of parameters. Name them
#'   (`prior_independent(beta = ..., gamma = ...)`) to name the parameters: a
#'   name given here beats the component's own `param_names` for a one-parameter
#'   component, and a multi-parameter component keeps its own names, or takes
#'   `name1`, `name2`, ... when it has none.
#' @return An `nsbi_prior` object.
#' @seealso [prior_families] for the components, [prior_truncated()] to bound
#'   one, [priors] for the rest.
#' @examples
#' prior <- prior_independent(
#'   p_S1S2 = prior_beta(2, 15),
#'   hr_S1  = prior_lognormal(log(3), 0.3),
#'   hr_S2  = prior_lognormal(log(10), 0.25)
#' )
#' prior
#' sample_prior(prior, 3)
#' @export
prior_independent <- function(...) {
  parts <- list(...)
  if (length(parts) == 0L) {
    stop("prior_independent() needs at least one prior. See ?priors.",
         call. = FALSE)
  }
  nms <- names(parts)
  for (i in seq_along(parts)) {
    if (inherits(parts[[i]], "nsbi_prior")) next
    where <- if (!is.null(nms) && nzchar(nms[[i]])) sprintf("`%s`", nms[[i]])
      else sprintf("argument %d", i)
    stop(sprintf(paste0("Every argument to prior_independent() must be an ",
                        "nsbi_prior object, but %s is of class %s. Build one ",
                        "with prior_uniform(), prior_normal() or one of the ",
                        "families in ?prior_families."),
                 where,
                 if (is.null(parts[[i]])) "NULL" else
                   paste(class(parts[[i]]), collapse = "/")),
         call. = FALSE)
  }

  dims <- vapply(parts, function(p) as.integer(p$dim), integer(1))
  d <- sum(dims)
  param_names <- independent_param_names(parts, dims, nms)

  # The fast path keeps the marginal specifications, and with them truncation
  # and the Stan export. A prior_custom() component has no specification, so
  # the whole product falls back to composing the closures.
  marginals <- lapply(parts, function(p) p$params$marginals)
  if (!any(vapply(marginals, is.null, logical(1)))) {
    return(marginal_prior(unlist(marginals, recursive = FALSE),
                          param_names = param_names, type = "independent"))
  }

  blocks <- split(seq_len(d), rep(seq_along(parts), dims))
  sample_fn <- function(n) {
    out <- do.call(cbind, lapply(parts, function(p) {
      as_theta_matrix(p$sample(n), p$dim)
    }))
    # Component column names would only ever be a partial set here; the joint
    # prior's own param_names are applied by sample_prior().
    dimnames(out) <- NULL
    out
  }
  log_prob_fn <- function(theta) {
    theta <- as_theta_matrix(theta, d)
    total <- rep(0, nrow(theta))
    for (i in seq_along(parts)) {
      total <- total +
        as.numeric(parts[[i]]$log_prob(theta[, blocks[[i]], drop = FALSE]))
    }
    total
  }
  lower <- unname(unlist(lapply(parts, function(p) p$lower %||% rep(-Inf, p$dim))))
  upper <- unname(unlist(lapply(parts, function(p) p$upper %||% rep(Inf, p$dim))))

  new_prior(sample_fn, log_prob_fn, d,
            lower = if (all(lower == -Inf)) NULL else lower,
            upper = if (all(upper == Inf)) NULL else upper,
            type = "independent", param_names = param_names,
            # The components themselves, not just their type names: this is
            # what lets is_improper_uniform_prior() (R/prior.R) recurse into
            # a composed prior instead of only checking the top-level type,
            # which is "independent" here regardless of what is inside.
            params = list(components = unname(parts)))
}

#' Parameter names for a product prior
#'
#' An argument name wins for a one-parameter component, since that is the
#' obvious way to write it; a wider component keeps the names it already has,
#' and otherwise the argument name is numbered. Anything short of a complete
#' set of names is dropped, because a partly named prior would label some
#' columns and not others.
#' @param parts The component priors.
#' @param dims Their dimensions.
#' @param nms The names of the `...` arguments, or `NULL`.
#' @keywords internal
independent_param_names <- function(parts, dims, nms) {
  out <- character(0)
  for (i in seq_along(parts)) {
    tag <- if (!is.null(nms) && nzchar(nms[[i]])) nms[[i]] else NA_character_
    own <- parts[[i]]$param_names
    got <- if (dims[[i]] == 1L) {
      if (!is.na(tag)) tag else own %||% NA_character_
    } else if (!is.null(own)) {
      own
    } else if (!is.na(tag)) {
      paste0(tag, seq_len(dims[[i]]))
    } else {
      rep(NA_character_, dims[[i]])
    }
    out <- c(out, got)
  }
  if (length(out) != sum(dims) || anyNA(out) || !all(nzchar(out))) NULL else out
}

#' Truncate a prior to a box
#'
#' Stan's `T[lower, upper]`, as an object. The returned prior is the original
#' restricted to the box and renormalized by the mass it keeps, so its
#' `log_prob` is a proper log density rather than the original shifted by an
#' unknown constant. That matters here in a way it does not in Stan: the
#' density is compared against a learned posterior in [nle()]'s MCMC target and
#' summed with it in [c2st()] and the diagnostics, so a missing constant is a
#' wrong answer rather than a constant offset.
#'
#' Renormalization needs the family's CDF, so `prior` has to come from a named
#' family: [prior_uniform()], [prior_normal()], anything in [prior_families],
#' or a [prior_independent()] built from those. A [prior_custom()] is arbitrary
#' R code with no CDF behind it; give it `lower`/`upper` there instead, which
#' rejects out-of-support draws without claiming to renormalize.
#'
#' @param prior An `nsbi_prior` from a named family.
#' @param lower,upper Truncation bounds. Numeric of length `prior$dim`, or
#'   length 1 to apply the same bound to every parameter. Give one or both;
#'   `-Inf`/`Inf` leaves that side alone.
#' @return An `nsbi_prior` object.
#' @seealso [prior_families], [prior_independent()], [priors].
#' @examples
#' # A half-normal, the long way round (prior_half_normal() is the short one).
#' prior_truncated(prior_normal(mean = 0, sd = 1), lower = 0)
#'
#' # A contact rate known to be between 0.1 and 2, log-normal in between.
#' prior <- prior_truncated(prior_lognormal(log(0.4), 0.5),
#'                          lower = 0.1, upper = 2)
#' range(sample_prior(prior, 100))
#' @export
prior_truncated <- function(prior, lower = NULL, upper = NULL) {
  check_prior(prior)
  d <- prior$dim
  marginals <- prior$params$marginals
  if (is.null(marginals)) {
    stop(sprintf(paste0("prior_truncated() needs a prior built from a named ",
                        "family, but this one has type \"%s\" and carries no ",
                        "family to renormalize against.\nUse the `lower`/",
                        "`upper` arguments of prior_custom() to bound it ",
                        "instead; they reject out-of-support draws without ",
                        "renormalizing the density. See ?prior_families."),
                 prior$type),
         call. = FALSE)
  }
  if (is.null(lower) && is.null(upper)) {
    stop("Give `lower`, `upper` or both; prior_truncated() has nothing to do ",
         "otherwise.", call. = FALSE)
  }
  lower <- check_bound(lower, "lower", d) %||% rep(-Inf, d)
  upper <- check_bound(upper, "upper", d) %||% rep(Inf, d)
  if (any(upper <= lower)) {
    stop("Every `upper` must be strictly greater than the matching `lower`.",
         call. = FALSE)
  }

  truncated <- lapply(seq_len(d), function(j) {
    m <- marginals[[j]]
    new_marginal(m$family, m$args,
                 lower = max(m$lower, lower[[j]]),
                 upper = min(m$upper, upper[[j]]),
                 label = param_label(prior$param_names, j))
  })
  marginal_prior(truncated, param_names = prior$param_names,
                 type = "truncated")
}

#' One line per marginal, for print.nsbi_prior()
#'
#' Reads back as close to the constructor call as one line allows:
#' `beta(2, 15) on [0, 1]`, with the bounds shown only when they are not the
#' family's own.
#' @param marginals The marginals to describe.
#' @keywords internal
describe_marginals <- function(marginals) {
  vapply(marginals, function(m) {
    args <- paste(vapply(m$args, function(v) format(signif(v, 4)),
                         character(1)), collapse = ", ")
    # Stan's own one-sided spelling: T[0, ] rather than T[0, Inf].
    trunc <- if (m$lower > m$natural[1L] || m$upper < m$natural[2L]) {
      sprintf(" T[%s, %s]",
              if (is.finite(m$lower)) format(m$lower) else "",
              if (is.finite(m$upper)) format(m$upper) else "")
    } else {
      ""
    }
    sprintf("%s(%s)%s", m$family, args, trunc)
  }, character(1))
}
