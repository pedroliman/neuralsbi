#' Neural Posterior Estimation (NPE)
#'
#' `npe()` is the main entry point. Given a prior and either a simulator (which
#' it will call) or a set of pre-computed simulations `(theta, x)`, it trains a
#' conditional density estimator whose output directly approximates the posterior
#' \eqn{p(\theta \mid x)}. This is single-round, *amortized* NPE: after training
#' once, you can condition on any observation without re-simulating.
#'
#' @param prior An `nsbi_prior` (see [prior_uniform()], [prior_normal()]).
#' @param simulator A function mapping an `n x dim` matrix of parameters to an
#'   `n x d` matrix of simulated data. Ignored if `theta` and `x` are given.
#' @param n_simulations Number of prior draws to simulate when `simulator` is
#'   used and `theta`/`x` are not supplied.
#' @param theta,x Optional pre-computed simulations. If supplied, `simulator`
#'   and `n_simulations` are ignored.
#' @param density_estimator One of `"mdn"` (neural Mixture Density Network,
#'   needs `torch`), `"maf"` (Masked Autoregressive Flow, needs `torch`),
#'   `"nsf"` (Neural Spline Flow, needs `torch`), or `"linear_gaussian"`
#'   (closed-form baseline, no `torch`), or a function `function(theta, x)`
#'   returning a fitted estimator.
#' @param n_transforms MAF/NSF setting: number of stacked autoregressive
#'   transforms.
#' @param n_components,hidden MDN settings: number of mixture components and a
#'   vector of hidden-layer widths.
#' @param max_epochs,batch_size,lr,validation_fraction,patience Neural training
#'   controls (Adam optimizer, early stopping on validation loss).
#' @param n_restarts Train this many independently initialized networks and keep
#'   the one with the best validation loss (guards against bad initializations
#'   and MDN mode collapse).
#' @param clip_grad_norm Maximum gradient norm during training (`Inf` disables
#'   clipping). The learning rate also decays 2x after 10 epochs without
#'   validation improvement.
#' @param standardize Whether to z-score `theta` and `x` before training
#'   (strongly recommended; default `TRUE`).
#' @param seed Optional integer seed for reproducibility.
#' @param verbose Print training progress.
#' @param ... Passed to the density estimator.
#'
#' @return An object of class `nsbi_npe`. Turn it into a usable posterior with
#'   [posterior()], or sample directly with [sample()].
#'
#' @examples
#' prior <- prior_uniform(c(-2, -2, -2), c(2, 2, 2))
#' simulator <- function(theta) theta + 1 + matrix(rnorm(length(theta), sd = 0.1),
#'                                                  nrow = nrow(theta))
#' fit <- npe(prior, simulator, n_simulations = 2000,
#'            density_estimator = "linear_gaussian")
#' post <- posterior(fit, x_obs = c(0.8, 0.6, 0.4))
#' draws <- sample(post, 1000)
#' @export
npe <- function(prior, simulator = NULL, n_simulations = 1000,
                theta = NULL, x = NULL,
                density_estimator = c("mdn", "maf", "nsf", "linear_gaussian"),
                n_components = 5L, n_transforms = 5L, hidden = c(50L, 50L),
                max_epochs = 500L, batch_size = 100L, lr = 5e-4,
                validation_fraction = 0.1, patience = 20L,
                n_restarts = 1L, clip_grad_norm = 5,
                standardize = TRUE, seed = NULL, verbose = FALSE, ...) {
  stopifnot(inherits(prior, "nsbi_prior"))

  if (is.null(theta) || is.null(x)) {
    if (is.null(simulator)) {
      stop("Provide either `simulator` or both `theta` and `x`.", call. = FALSE)
    }
    sims <- simulate_for_sbi(simulator, prior, n_simulations, seed = seed,
                             verbose = verbose)
    theta <- sims$theta
    x <- sims$x
  }
  theta <- as_theta_matrix(theta, prior$dim)
  x <- as_theta_matrix(x)
  if (nrow(theta) != nrow(x)) {
    stop("`theta` and `x` must have the same number of rows.", call. = FALSE)
  }

  # standardization
  if (standardize) {
    std_theta <- fit_standardizer(theta)
    std_x <- fit_standardizer(x)
  } else {
    std_theta <- fit_standardizer(matrix(0, 1, ncol(theta)))
    std_theta$scale[] <- 1; std_theta$center[] <- 0
    std_x <- fit_standardizer(matrix(0, 1, ncol(x)))
    std_x$scale[] <- 1; std_x$center[] <- 0
  }
  theta_z <- apply_standardizer(std_theta, theta)
  x_z <- apply_standardizer(std_x, x)

  de <- fit_density_estimator(
    density_estimator, theta_z, x_z,
    n_components = n_components, n_transforms = n_transforms,
    hidden = hidden, max_epochs = max_epochs,
    batch_size = batch_size, lr = lr, validation_fraction = validation_fraction,
    patience = patience, n_restarts = n_restarts,
    clip_grad_norm = clip_grad_norm, seed = seed, verbose = verbose, ...
  )

  structure(
    list(
      de = de,
      prior = prior,
      std_theta = std_theta,
      std_x = std_x,
      dim_theta = prior$dim,
      dim_x = ncol(x),
      n_simulations = nrow(theta),
      density_estimator = if (is.character(density_estimator))
        density_estimator[1] else "custom"
    ),
    class = "nsbi_npe"
  )
}

#' @keywords internal
fit_density_estimator <- function(density_estimator, theta_z, x_z, ...) {
  if (is.function(density_estimator)) {
    return(density_estimator(theta_z, x_z))
  }
  density_estimator <- match.arg(density_estimator,
                                 c("mdn", "maf", "nsf", "linear_gaussian"))
  dots <- list(...)
  switch(
    density_estimator,
    mdn = fit_mdn(theta_z, x_z,
                  n_components = dots$n_components %||% 5L,
                  hidden = dots$hidden %||% c(50L, 50L),
                  max_epochs = dots$max_epochs %||% 500L,
                  batch_size = dots$batch_size %||% 100L,
                  lr = dots$lr %||% 5e-4,
                  validation_fraction = dots$validation_fraction %||% 0.1,
                  patience = dots$patience %||% 20L,
                  n_restarts = dots$n_restarts %||% 1L,
                  clip_grad_norm = dots$clip_grad_norm %||% 5,
                  seed = dots$seed, verbose = dots$verbose %||% FALSE),
    maf = fit_maf(theta_z, x_z,
                  n_transforms = dots$n_transforms %||% 5L,
                  hidden = dots$hidden %||% c(50L, 50L),
                  max_epochs = dots$max_epochs %||% 500L,
                  batch_size = dots$batch_size %||% 100L,
                  lr = dots$lr %||% 5e-4,
                  validation_fraction = dots$validation_fraction %||% 0.1,
                  patience = dots$patience %||% 20L,
                  n_restarts = dots$n_restarts %||% 1L,
                  clip_grad_norm = dots$clip_grad_norm %||% 5,
                  seed = dots$seed, verbose = dots$verbose %||% FALSE),
    nsf = fit_nsf(theta_z, x_z,
                  n_transforms = dots$n_transforms %||% 5L,
                  hidden = dots$hidden %||% c(50L, 50L),
                  n_bins = dots$n_bins %||% 8L,
                  tail_bound = dots$tail_bound %||% 3,
                  max_epochs = dots$max_epochs %||% 500L,
                  batch_size = dots$batch_size %||% 100L,
                  lr = dots$lr %||% 5e-4,
                  validation_fraction = dots$validation_fraction %||% 0.1,
                  patience = dots$patience %||% 20L,
                  n_restarts = dots$n_restarts %||% 1L,
                  clip_grad_norm = dots$clip_grad_norm %||% 5,
                  seed = dots$seed, verbose = dots$verbose %||% FALSE),
    linear_gaussian = fit_linear_gaussian(theta_z, x_z,
                                           verbose = dots$verbose %||% FALSE)
  )
}

#' Run a simulator over prior draws
#'
#' @inheritParams npe
#' @param n Number of simulations.
#' @return A list with `theta` (`n x dim`) and `x` (`n x d`) matrices.
#' @export
simulate_for_sbi <- function(simulator, prior, n, seed = NULL, verbose = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  theta <- sample_prior(prior, n)
  verbose_cat(verbose, sprintf("Simulating %d draws...\n", n))
  x <- simulator(theta)
  x <- as_theta_matrix(x)
  if (nrow(x) != nrow(theta)) {
    stop("Simulator must return one row of output per row of theta.",
         call. = FALSE)
  }
  list(theta = theta, x = x)
}

#' @export
print.nsbi_npe <- function(x, ...) {
  cat("<nsbi_npe> Neural Posterior Estimation fit\n")
  cat(sprintf("  density estimator : %s\n", x$density_estimator))
  cat(sprintf("  parameters (dim)  : %d\n", x$dim_theta))
  cat(sprintf("  data (dim)        : %d\n", x$dim_x))
  cat(sprintf("  simulations       : %d\n", x$n_simulations))
  if (!is.null(x$de$best_val_loss) && is.finite(x$de$best_val_loss)) {
    cat(sprintf("  best val loss     : %.4f\n", x$de$best_val_loss))
  }
  cat("  -> build a posterior with posterior(fit, x_obs = ...)\n")
  invisible(x)
}
