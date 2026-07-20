#' @import S7
NULL

# ===========================================================================
# S7 class + generic definitions for neuralsbi.
#
# This file is named `aaa-classes.R` so it is collated first: every S7 class
# and generic is defined here, and the `method(<generic>, <class>) <- ...`
# registrations live in the topical files (prior.R, mdn.R, ...), which are
# collated afterwards. Keeping all type definitions in one place makes the
# object model easy to read and guarantees classes/generics exist before any
# method is registered.
# ===========================================================================

# ---- Priors ---------------------------------------------------------------

#' Prior distribution classes
#'
#' `Prior` is the abstract base class for priors. Concrete subclasses are
#' [UniformPrior], [NormalPrior] and [CustomPrior]; construct them with the
#' user-facing helpers [prior_uniform()], [prior_normal()], [prior_custom()].
#' Every prior carries its dimensionality and (optionally) box support bounds.
#'
#' @param n_dim Integer number of parameters.
#' @param lower,upper Optional numeric support bounds (or `NULL` if unbounded).
#' @export
Prior <- new_class(
  "Prior", abstract = TRUE,
  properties = list(
    n_dim = class_integer,
    # class_any (not `class_numeric | NULL`): an S7 union default of NULL
    # materializes as numeric(0), which would make unbounded priors look bounded.
    lower = new_property(class_any, default = NULL),
    upper = new_property(class_any, default = NULL)
  )
)

#' @rdname Prior
#' @param low,high Numeric vectors of box bounds.
#' @export
UniformPrior <- new_class(
  "UniformPrior", parent = Prior,
  properties = list(low = class_numeric, high = class_numeric)
)

#' @rdname Prior
#' @param mean Numeric vector of means.
#' @param sd Numeric vector of standard deviations.
#' @export
NormalPrior <- new_class(
  "NormalPrior", parent = Prior,
  properties = list(mean = class_numeric, sd = class_numeric)
)

#' @rdname Prior
#' @param sample_fn Function `function(n)` returning an `n x n_dim` matrix.
#' @param log_prob_fn Function `function(theta)` returning a length-`n` vector.
#' @export
CustomPrior <- new_class(
  "CustomPrior", parent = Prior,
  properties = list(sample_fn = class_function, log_prob_fn = class_function)
)

#' Prior generics
#'
#' * `draw_prior(prior, n)` -- draw `n` samples (`n x n_dim` matrix).
#' * `prior_log_prob(prior, theta)` -- log density, one value per row.
#' * `in_support(prior, theta)` -- logical, whether each row is in support.
#'
#' @param prior A [Prior] object.
#' @param n Number of draws.
#' @param theta Matrix (or vector) of parameters.
#' @name prior-generics
#' @export
draw_prior <- new_generic("draw_prior", "prior")
#' @rdname prior-generics
#' @export
prior_log_prob <- new_generic("prior_log_prob", "prior")
#' @rdname prior-generics
#' @export
in_support <- new_generic("in_support", "prior")

# ---- Density estimators ---------------------------------------------------

#' Conditional density estimator classes
#'
#' `DensityEstimator` is the abstract base for estimators of \eqn{q(\theta\mid x)}.
#' [MDN] is a neural Mixture Density Network (needs `torch`); [LinearGaussian] is
#' a closed-form conditional Gaussian baseline. Both train in standardized space
#' and implement the generics [de_log_prob()] and [de_sample()].
#'
#' @export
DensityEstimator <- new_class(
  "DensityEstimator", abstract = TRUE,
  properties = list(n_dim_theta = class_integer, n_dim_x = class_integer)
)

#' @rdname DensityEstimator
#' @export
MDN <- new_class(
  "MDN", parent = DensityEstimator,
  properties = list(
    net = class_any,
    n_components = class_integer,
    hidden = class_numeric,
    best_val_loss = new_property(class_numeric, default = NA_real_)
  )
)

#' @rdname DensityEstimator
#' @export
LinearGaussian <- new_class(
  "LinearGaussian", parent = DensityEstimator,
  properties = list(B = class_any, Sigma = class_any, chol = class_any)
)

#' Density-estimator generics
#'
#' @param de A [DensityEstimator].
#' @param theta,x Matrices of parameters / data (standardized space).
#' @param n Number of draws.
#' @name de-generics
#' @export
de_log_prob <- new_generic("de_log_prob", "de")
#' @rdname de-generics
#' @export
de_sample <- new_generic("de_sample", "de")

# ---- Trainer (sbi-style builder) ------------------------------------------

#' Neural Posterior Estimation trainer (S7)
#'
#' `NPE` mirrors the builder workflow of the Python `sbi` package:
#' `NPE(prior)` -> [append_simulations()] -> [train()] -> [build_posterior()].
#' Each step returns an updated object. For a one-call interface use [npe()].
#'
#' @param prior A [Prior].
#' @param density_estimator Estimator name, `"mdn"` or `"linear_gaussian"`.
#' @param config List of estimator/training hyperparameters.
#' @export
NPE <- new_class(
  "NPE",
  properties = list(
    prior = Prior,
    theta = new_property(class_any, default = NULL),
    x = new_property(class_any, default = NULL),
    # class_any (not `DensityEstimator | NULL`): an S7 union over an abstract
    # class tries to construct the abstract prototype. Checked in methods.
    estimator = new_property(class_any, default = NULL),
    std_theta = new_property(class_any, default = NULL),
    std_x = new_property(class_any, default = NULL),
    density_estimator = new_property(class_character, default = "mdn"),
    config = new_property(class_list, default = list())
  )
)

#' NPE builder generics
#'
#' @param npe An [NPE] object.
#' @param theta,x Simulated parameters and data.
#' @param x_obs Optional default observation for the posterior.
#' @param ... Passed on.
#' @name npe-generics
#' @export
append_simulations <- new_generic("append_simulations", "npe")
#' @rdname npe-generics
#' @export
train <- new_generic("train", "npe")
#' @rdname npe-generics
#' @export
build_posterior <- new_generic("build_posterior", "npe")

# ---- Posteriors -----------------------------------------------------------

#' Posterior classes
#'
#' `Posterior` is the abstract base. NPE yields a [DirectPosterior]: the trained
#' density estimator *is* the posterior, so sampling and density evaluation are
#' immediate (no MCMC). It wraps the [NPE] fit plus an optional default
#' observation.
#'
#' @export
Posterior <- new_class("Posterior", abstract = TRUE)

#' @rdname Posterior
#' @export
DirectPosterior <- new_class(
  "DirectPosterior", parent = Posterior,
  properties = list(
    fit = NPE,
    default_x = new_property(class_any, default = NULL)
  )
)

#' Posterior generics
#'
#' @param x,post A [Posterior] object (`x` for the [sample()] generic).
#' @param theta Parameters to evaluate.
#' @param obs,x Observation to condition on.
#' @param ... Passed on.
#' @name posterior-generics
#' @export
sample <- new_generic("sample", "x")
#' @rdname posterior-generics
#' @export
log_prob <- new_generic("log_prob", "post")
#' @rdname posterior-generics
#' @export
map_estimate <- new_generic("map_estimate", "post")

# ---- Diagnostics result ---------------------------------------------------

#' Simulation-based calibration result
#'
#' Value object returned by [sbc()]; consumed by [expected_coverage()] and
#' [plot_sbc()].
#'
#' @export
SBCResult <- new_class(
  "SBCResult",
  properties = list(
    ranks = class_any,
    n_posterior_samples = class_integer,
    n_sbc = class_integer,
    uniformity_pvalue = class_numeric
  )
)
