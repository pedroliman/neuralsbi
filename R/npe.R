#' Neural Posterior Estimation (NPE)
#'
#' `npe()` is the one-call entry point. Given a prior and either a simulator
#' (which it will run) or pre-computed simulations `(theta, x)`, it trains a
#' conditional density estimator whose output directly approximates the
#' posterior \eqn{p(\theta\mid x)}. This is single-round, *amortized* NPE: train
#' once, then condition on any observation without re-simulating.
#'
#' Internally this builds an [NPE] object and runs
#' [append_simulations()] |> [train()]; use those directly for the sbi-style
#' builder workflow.
#'
#' @param prior A [Prior] (see [prior_uniform()], [prior_normal()]).
#' @param simulator A function mapping an `n x n_dim` parameter matrix to an
#'   `n x d` data matrix. Ignored if `theta` and `x` are given.
#' @param n_simulations Number of prior draws to simulate when `simulator` is
#'   used and `theta`/`x` are not supplied.
#' @param theta,x Optional pre-computed simulations.
#' @param density_estimator `"mdn"` (neural, needs `torch`) or
#'   `"linear_gaussian"` (closed-form baseline), or a function
#'   `function(theta, x)` returning a fitted [DensityEstimator].
#' @param n_components,hidden MDN mixture components and hidden-layer widths.
#' @param max_epochs,batch_size,lr,validation_fraction,patience MDN training
#'   controls.
#' @param standardize Whether to z-score `theta` and `x` (recommended).
#' @param seed Optional integer seed.
#' @param verbose Print progress.
#' @param ... Passed to the density estimator.
#' @return A trained [NPE] object. Build a posterior with [posterior()].
#' @examples
#' \donttest{
#' prior <- prior_uniform(c(-2, -2, -2), c(2, 2, 2))
#' simulator <- function(theta) theta + 1 +
#'   matrix(rnorm(length(theta), sd = 0.1), nrow = nrow(theta))
#' fit <- npe(prior, simulator, n_simulations = 2000)
#' post <- posterior(fit, x_obs = c(0.8, 0.6, 0.4))
#' draws <- sample(post, 1000)
#' }
#' @export
npe <- function(prior, simulator = NULL, n_simulations = 1000,
                theta = NULL, x = NULL,
                density_estimator = c("mdn", "linear_gaussian"),
                n_components = 5L, hidden = c(50L, 50L),
                max_epochs = 500L, batch_size = 100L, lr = 5e-4,
                validation_fraction = 0.1, patience = 20L,
                standardize = TRUE, seed = NULL, verbose = FALSE, ...) {
  de_spec <- if (is.function(density_estimator)) density_estimator else
    match.arg(density_estimator)
  if (is.null(theta) || is.null(x)) {
    if (is.null(simulator)) {
      stop("Provide either `simulator` or both `theta` and `x`.", call. = FALSE)
    }
    sims <- simulate_for_sbi(simulator, prior, n_simulations, seed = seed,
                             verbose = verbose)
    theta <- sims$theta; x <- sims$x
  }
  trainer <- NPE(
    prior = prior,
    density_estimator = if (is.function(de_spec)) "custom" else de_spec,
    config = list(density_estimator = de_spec, n_components = n_components,
                  hidden = hidden, max_epochs = max_epochs,
                  batch_size = batch_size, lr = lr,
                  validation_fraction = validation_fraction,
                  patience = patience, standardize = standardize,
                  seed = seed, verbose = verbose, dots = list(...))
  )
  trainer <- append_simulations(trainer, theta, x)
  train(trainer)
}

# ---- builder methods (sbi-style) ------------------------------------------

method(append_simulations, NPE) <- function(npe, theta, x) {
  theta <- as_theta_matrix(theta, npe@prior@n_dim)
  x <- as_theta_matrix(x)
  if (nrow(theta) != nrow(x)) {
    stop("`theta` and `x` must have the same number of rows.", call. = FALSE)
  }
  if (!is.null(npe@theta)) {
    theta <- rbind(npe@theta, theta)
    x <- rbind(npe@x, x)
  }
  npe@theta <- theta
  npe@x <- x
  npe
}

method(train, NPE) <- function(npe, ...) {
  cfg <- utils::modifyList(npe@config, list(...))
  theta <- npe@theta; x <- npe@x
  if (is.null(theta) || is.null(x)) {
    stop("No simulations. Call append_simulations() first.", call. = FALSE)
  }
  standardize <- cfg$standardize %||% TRUE
  if (standardize) {
    std_theta <- fit_standardizer(theta)
    std_x <- fit_standardizer(x)
  } else {
    std_theta <- identity_standardizer(ncol(theta))
    std_x <- identity_standardizer(ncol(x))
  }
  theta_z <- apply_standardizer(std_theta, theta)
  x_z <- apply_standardizer(std_x, x)

  # estimator spec: from config (functional npe() path) or the @density_estimator
  # slot (sbi-style builder path).
  de_spec <- cfg$density_estimator %||% npe@density_estimator
  de <- fit_density_estimator(de_spec, theta_z, x_z, cfg)

  npe@estimator <- de
  npe@std_theta <- std_theta
  npe@std_x <- std_x
  npe
}

method(build_posterior, NPE) <- function(npe, x_obs = NULL, ...) {
  if (is.null(npe@estimator)) stop("Train the NPE first.", call. = FALSE)
  if (!is.null(x_obs)) x_obs <- as_theta_matrix(x_obs, ncol(npe@x))
  DirectPosterior(fit = npe, default_x = x_obs)
}

#' @keywords internal
fit_density_estimator <- function(density_estimator, theta_z, x_z, cfg) {
  if (is.function(density_estimator)) return(density_estimator(theta_z, x_z))
  switch(
    density_estimator,
    mdn = fit_mdn(theta_z, x_z,
                  n_components = cfg$n_components %||% 5L,
                  hidden = cfg$hidden %||% c(50L, 50L),
                  max_epochs = cfg$max_epochs %||% 500L,
                  batch_size = cfg$batch_size %||% 100L,
                  lr = cfg$lr %||% 5e-4,
                  validation_fraction = cfg$validation_fraction %||% 0.1,
                  patience = cfg$patience %||% 20L,
                  seed = cfg$seed, verbose = cfg$verbose %||% FALSE),
    linear_gaussian = fit_linear_gaussian(theta_z, x_z,
                                          verbose = cfg$verbose %||% FALSE),
    stop("Unknown density_estimator '", density_estimator, "'.", call. = FALSE)
  )
}

#' Run a simulator over prior draws
#'
#' @inheritParams npe
#' @param n Number of simulations.
#' @return A list with `theta` (`n x n_dim`) and `x` (`n x d`) matrices.
#' @export
simulate_for_sbi <- function(simulator, prior, n, seed = NULL, verbose = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  theta <- sample_prior(prior, n)
  verbose_cat(verbose, sprintf("Simulating %d draws...\n", n))
  x <- as_theta_matrix(simulator(theta))
  if (nrow(x) != nrow(theta)) {
    stop("Simulator must return one row of output per row of theta.",
         call. = FALSE)
  }
  list(theta = theta, x = x)
}

#' Number of simulations a trained NPE has seen
#' @param fit An [NPE] object.
#' @return Integer count of simulations.
#' @export
n_simulations <- function(fit) {
  if (is.null(fit@theta)) 0L else nrow(fit@theta)
}

method(print, NPE) <- function(x, ...) {
  cat("<NPE> Neural Posterior Estimation\n")
  cat(sprintf("  density estimator : %s\n", x@density_estimator))
  cat(sprintf("  parameters (dim)  : %d\n", x@prior@n_dim))
  cat(sprintf("  simulations       : %d\n", n_simulations(x)))
  if (is.null(x@estimator)) {
    cat("  status            : not trained (call train())\n")
  } else {
    cat(sprintf("  data (dim)        : %d\n", x@estimator@n_dim_x))
    if (S7_inherits(x@estimator, MDN) && is.finite(x@estimator@best_val_loss)) {
      cat(sprintf("  best val loss     : %.4f\n", x@estimator@best_val_loss))
    }
    cat("  -> posterior(fit, x_obs = ...)\n")
  }
  invisible(x)
}
