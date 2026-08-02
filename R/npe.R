#' Neural Posterior Estimation (NPE)
#'
#' `npe()` is the main entry point. Given a prior and either a simulator (which
#' it will call) or a set of pre-computed simulations `(theta, x)`, it trains a
#' conditional density estimator whose output directly approximates the posterior
#' \eqn{p(\theta \mid x)}. This is single-round, *amortized* NPE: after training
#' once, you can condition on any observation without re-simulating.
#'
#' @param prior An `nsbi_prior` (see [prior_uniform()], [prior_normal()]).
#' @param simulator A function called once per parameter set, returning one
#'   simulated observation: a numeric vector, a scalar, or a one-row matrix or
#'   data frame. Parameters arrive either as named arguments (when the prior's
#'   names match the simulator's formals) or as one named vector. Names on the
#'   output become the outcome names used in plots. See [nsbi_simulator].
#'   Ignored if `theta` and `x` are given.
#' @param n_simulations Number of prior draws to simulate when `simulator` is
#'   used and `theta`/`x` are not supplied.
#' @param sim_args Named list of extra arguments passed to every simulator
#'   call: observed data, a time grid, a fixed population size, solver
#'   settings. See [nsbi_simulator].
#' @param theta,x Optional pre-computed simulations. If supplied, `simulator`
#'   and `n_simulations` are ignored. Column names on `theta` (or names on
#'   `prior`'s `mean`/`low`) and on `x` are carried through to posterior
#'   samples, SBC results, and their plots.
#' @param density_estimator One of `"maf"` (Masked Autoregressive Flow, needs
#'   `torch`; the default, matching Python `sbi`), `"mdn"` (neural Mixture
#'   Density Network, needs `torch`), `"nsf"` (Neural Spline Flow, needs
#'   `torch`), or `"linear_gaussian"` (closed-form baseline, no `torch`), or a
#'   function `function(theta, x)` returning a fitted estimator.
#' @param n_transforms MAF/NSF setting: number of stacked autoregressive
#'   transforms (default 5, as in `sbi`).
#' @param n_components,hidden MDN settings: number of mixture components
#'   (default 10, as in `sbi`) and a vector of hidden-layer widths.
#' @param n_bins,tail_bound NSF settings: number of spline bins per transform
#'   and the half-width of the interval the spline acts on (outside it the
#'   transform is the identity).
#' @param embedding_net Optional summary network built with [embedding_mlp()].
#'   When supplied, the neural estimators condition on the learned features
#'   \eqn{f_\psi(x)} instead of the raw data, training the embedding jointly.
#'   Ignored (with a warning) by `"linear_gaussian"`.
#' @param max_epochs,batch_size,lr,validation_fraction,patience Neural training
#'   controls (Adam optimizer, early stopping on validation loss). The defaults
#'   (`batch_size = 200`, `lr = 5e-4`, `validation_fraction = 0.1`,
#'   `patience = 20`) match Python `sbi`; `max_epochs` is a high guard cap that
#'   early stopping normally reaches first.
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
#'
#' @section Parallel simulation and progress:
#'
#' The simulator runs sequentially unless you declare a \pkg{future} plan --
#' `library(future); plan(multisession)` -- in which case the simulations are
#' spread across workers. Simulation and training both report progress
#' with an ETA. See [nsbi_parallel] and [nsbi_progress].
#'
#' @return An object of class `nsbi_npe`. Turn it into a usable posterior with
#'   [posterior()], or sample directly with [sample()]. Save it to disk with
#'   [save_npe()]: a torch-backed fit does not survive `saveRDS()`.
#'
#' @examples
#' prior <- prior_uniform(c(mu = -2, nu = -2), c(mu = 2, nu = 2))
#' simulator <- function(mu, nu) c(a = mu + rnorm(1, sd = 0.1),
#'                                 b = nu + rnorm(1, sd = 0.1))
#' fit <- npe(prior, simulator, n_simulations = 2000,
#'            density_estimator = "linear_gaussian")
#' post <- posterior(fit, x_obs = c(0.8, 0.6))
#' draws <- sample(post, 1000)
#' @export
npe <- function(prior, simulator = NULL, n_simulations = 1000,
                sim_args = list(), theta = NULL, x = NULL,
                density_estimator = c("maf", "mdn", "nsf", "linear_gaussian"),
                n_components = 10L, n_transforms = 5L, hidden = c(50L, 50L),
                n_bins = 10L, tail_bound = 3,
                embedding_net = NULL,
                max_epochs = 2000L, batch_size = 200L, lr = 5e-4,
                validation_fraction = 0.1, patience = 20L,
                n_restarts = 1L, clip_grad_norm = 5,
                standardize = TRUE, seed = NULL, verbose = FALSE) {
  # Everything here runs before the simulator does. An argument that is only
  # noticed by the arithmetic inside training costs the whole budget first,
  # and the budget is the expensive part of a run.
  if (!is.function(density_estimator)) {
    density_estimator <- match.arg(density_estimator)
  }
  check_prior(prior)
  if (!is.null(embedding_net) && !inherits(embedding_net, "nsbi_embedding")) {
    stop("`embedding_net` must be built with embedding_mlp().", call. = FALSE)
  }
  # tested against the matched name, so `density_estimator = "linear"` warns too
  if (!is.null(embedding_net) && identical(density_estimator,
                                           "linear_gaussian")) {
    warning("`embedding_net` is ignored by the linear_gaussian estimator.",
            call. = FALSE)
  }
  check_architecture(n_components, n_transforms, hidden, n_bins, tail_bound)
  check_train_controls(max_epochs, batch_size, lr, validation_fraction,
                       patience, n_restarts, clip_grad_norm)

  prep <- prepare_simulations(prior, simulator, n_simulations, sim_args,
                              theta, x, standardize, seed, verbose)
  theta <- prep$theta
  x <- prep$x
  n_dropped <- prep$n_dropped
  param_names <- prep$param_names
  x_names <- prep$x_names
  std_theta <- prep$std_theta
  std_x <- prep$std_x
  theta_z <- prep$theta_z
  x_z <- prep$x_z

  de <- fit_density_estimator(
    density_estimator, theta_z, x_z,
    n_components = n_components, n_transforms = n_transforms,
    hidden = hidden, n_bins = n_bins, tail_bound = tail_bound,
    embedding_net = embedding_net, max_epochs = max_epochs,
    batch_size = batch_size, lr = lr, validation_fraction = validation_fraction,
    patience = patience, n_restarts = n_restarts,
    clip_grad_norm = clip_grad_norm, seed = seed, verbose = verbose
  )

  structure(
    list(
      de = de,
      prior = prior,
      std_theta = std_theta,
      std_x = std_x,
      dim_theta = prior$dim,
      dim_x = ncol(x),
      param_names = param_names,
      x_names = x_names,
      n_simulations = nrow(theta),
      n_dropped = n_dropped,
      density_estimator = if (is.character(density_estimator))
        density_estimator[1] else "custom"
    ),
    class = "nsbi_npe"
  )
}

#' Everything both [npe()] and [nle()] do before touching a density estimator
#'
#' Get the simulations (running the simulator or taking pre-computed ones),
#' coerce them to matrices, drop non-finite draws, and learn the two
#' standardizers. Shared so the two entry points cannot drift apart -- they
#' differ only in which side of `(theta, x)` the estimator treats as its target.
#'
#' @return A list with the cleaned `theta`/`x`, their standardized versions,
#'   the two `nsbi_standardizer`s, the column names, and `n_dropped`.
#' @keywords internal
prepare_simulations <- function(prior, simulator, n_simulations, sim_args,
                                theta, x, standardize, seed, verbose) {
  if (is.null(theta) || is.null(x)) {
    if (is.null(simulator)) {
      stop("Provide either `simulator` or both `theta` and `x`.", call. = FALSE)
    }
    n_simulations <- check_count(
      n_simulations, "n_simulations", min = 2L,
      why = paste("since the estimator is standardized and split into a",
                  "training and a validation set"))
    sims <- simulate_for_sbi(simulator, prior, n_simulations,
                             sim_args = sim_args, seed = seed,
                             verbose = verbose)
    theta <- sims$theta
    x <- sims$x
    n_dropped <- sims$n_dropped
  } else {
    # pre-computed simulations get the same type and finiteness checks, so the
    # rules do not depend on who ran the simulator
    theta <- as_theta_matrix(check_numeric(theta, "theta"), prior$dim)
    x <- as_theta_matrix(check_numeric(x, "x"))
    if (nrow(theta) != nrow(x)) {
      stop("`theta` and `x` must have the same number of rows.", call. = FALSE)
    }
    kept <- drop_failed_sims(theta, x,
                             x_hint = "Check `x` for NA, NaN or Inf.")
    theta <- kept$theta
    x <- kept$x
    n_dropped <- kept$n_dropped
  }
  theta <- as_theta_matrix(theta, prior$dim)
  x <- as_theta_matrix(x)

  if (standardize) {
    std_theta <- fit_standardizer(theta, what = "theta")
    std_x <- fit_standardizer(x, what = "x")
  } else {
    std_theta <- fit_standardizer(matrix(0, 1, ncol(theta)))
    std_theta$scale[] <- 1; std_theta$center[] <- 0
    std_x <- fit_standardizer(matrix(0, 1, ncol(x)))
    std_x$scale[] <- 1; std_x$center[] <- 0
  }

  list(
    theta = theta, x = x,
    theta_z = apply_standardizer(std_theta, theta),
    x_z = apply_standardizer(std_x, x),
    std_theta = std_theta, std_x = std_x,
    param_names = colnames(theta) %||% prior$param_names,
    x_names = colnames(x),
    n_dropped = n_dropped
  )
}

#' Validate the architecture arguments shared by [npe()] and [nle()]
#'
#' All of them are checked whichever estimator was asked for, so `n_bins = 0`
#' is an error even under `"linear_gaussian"`, which ignores it. A value that
#' cannot build a network is a mistake in the call whether or not this run
#' would have read it.
#'
#' @inheritParams npe
#' @keywords internal
check_architecture <- function(n_components, n_transforms, hidden, n_bins,
                               tail_bound) {
  check_count(n_components, "n_components")
  check_count(n_transforms, "n_transforms")
  check_counts(hidden, "hidden", what = "one hidden-layer width per entry")
  check_count(n_bins, "n_bins")
  check_positive(tail_bound, "tail_bound")
  invisible(TRUE)
}

#' @keywords internal
fit_density_estimator <- function(density_estimator, theta_z, x_z, ...) {
  if (is.function(density_estimator)) {
    return(density_estimator(theta_z, x_z))
  }
  density_estimator <- match.arg(density_estimator,
                                 c("maf", "mdn", "nsf", "linear_gaussian"))
  dots <- list(...)
  switch(
    density_estimator,
    mdn = fit_mdn(theta_z, x_z,
                  n_components = dots$n_components %||% 10L,
                  hidden = dots$hidden %||% c(50L, 50L),
                  max_epochs = dots$max_epochs %||% 2000L,
                  batch_size = dots$batch_size %||% 200L,
                  lr = dots$lr %||% 5e-4,
                  validation_fraction = dots$validation_fraction %||% 0.1,
                  patience = dots$patience %||% 20L,
                  n_restarts = dots$n_restarts %||% 1L,
                  clip_grad_norm = dots$clip_grad_norm %||% 5,
                  embedding = dots$embedding_net,
                  seed = dots$seed, verbose = dots$verbose %||% FALSE),
    maf = fit_maf(theta_z, x_z,
                  n_transforms = dots$n_transforms %||% 5L,
                  hidden = dots$hidden %||% c(50L, 50L),
                  max_epochs = dots$max_epochs %||% 2000L,
                  batch_size = dots$batch_size %||% 200L,
                  lr = dots$lr %||% 5e-4,
                  validation_fraction = dots$validation_fraction %||% 0.1,
                  patience = dots$patience %||% 20L,
                  n_restarts = dots$n_restarts %||% 1L,
                  clip_grad_norm = dots$clip_grad_norm %||% 5,
                  embedding = dots$embedding_net,
                  seed = dots$seed, verbose = dots$verbose %||% FALSE),
    nsf = fit_nsf(theta_z, x_z,
                  n_transforms = dots$n_transforms %||% 5L,
                  hidden = dots$hidden %||% c(50L, 50L),
                  n_bins = dots$n_bins %||% 10L,
                  tail_bound = dots$tail_bound %||% 3,
                  max_epochs = dots$max_epochs %||% 2000L,
                  batch_size = dots$batch_size %||% 200L,
                  lr = dots$lr %||% 5e-4,
                  validation_fraction = dots$validation_fraction %||% 0.1,
                  patience = dots$patience %||% 20L,
                  n_restarts = dots$n_restarts %||% 1L,
                  clip_grad_norm = dots$clip_grad_norm %||% 5,
                  embedding = dots$embedding_net,
                  seed = dots$seed, verbose = dots$verbose %||% FALSE),
    linear_gaussian = fit_linear_gaussian(theta_z, x_z,
                                           verbose = dots$verbose %||% FALSE)
  )
}

#' Run a simulator over prior draws
#'
#' Draws `n` parameter vectors from the prior and calls the simulator once per
#' draw. Under a \pkg{future} plan the draws are spread across workers. See
#' [nsbi_simulator], [nsbi_parallel] and [nsbi_progress].
#'
#' Simulations whose output is not finite are dropped together with their
#' parameters, with a warning.
#'
#' The simulator comes first here and second in [npe()], [nle()] and
#' [npe_sequential()]. That is the reverse of the fitting functions and it is
#' easy to get backwards, so a call with the two swapped is detected and named
#' rather than left to fail inside [sample_prior()].
#'
#' @inheritParams npe
#' @param simulator A function called once per parameter set, returning one
#'   simulated observation: a numeric vector, a scalar, or a one-row matrix or
#'   data frame. See [nsbi_simulator]. Note the order: the simulator is the
#'   first argument here and the second in `npe(prior, simulator, ...)`.
#' @param n Number of simulations.
#' @return A list with `theta` (`n x dim`) and `x` (`n x d`) matrices and
#'   `n_dropped`, the number of simulations discarded for non-finite output.
#' @examples
#' prior <- prior_uniform(c(a = -1, b = -1), c(a = 1, b = 1))
#' sims <- simulate_for_sbi(function(a, b) c(a^2, b^2), prior, n = 100)
#' str(sims)
#' @export
simulate_for_sbi <- function(simulator, prior, n, sim_args = list(),
                             seed = NULL, verbose = FALSE) {
  # A prior object is never a function and a simulator is never an nsbi_prior,
  # so this pattern can only be the arguments the wrong way round. It goes
  # first: the general checks below would blame whichever argument they reached
  # and send the user looking at a prior that is fine.
  if (inherits(simulator, "nsbi_prior") && is.function(prior)) {
    stop("`simulator` and `prior` look swapped. simulate_for_sbi() takes the ",
         "simulator first: simulate_for_sbi(simulator, prior, n). Note this ",
         "is the reverse of npe(prior, simulator, ...).", call. = FALSE)
  }
  check_function(simulator, "simulator", what = "one parameter set per call")
  check_prior(prior)
  n <- check_count(n, "n")
  if (!is.null(seed)) set.seed(seed)
  theta <- sample_prior(prior, n)
  verbose_cat(verbose, sprintf("Simulating %d draws...\n", n))
  x <- run_simulator(simulator, theta, sim_args = sim_args)
  drop_failed_sims(theta, x)[c("theta", "x", "n_dropped")]
}

#' @export
print.nsbi_npe <- function(x, ...) {
  cat("<nsbi_npe> Neural Posterior Estimation fit\n")
  cat(sprintf("  density estimator : %s\n", x$density_estimator))
  cat(sprintf("  parameters (dim)  : %d\n", x$dim_theta))
  if (!is.null(x$param_names)) {
    cat("    names           :", paste(x$param_names, collapse = ", "), "\n")
  }
  cat(sprintf("  data (dim)        : %d\n", x$dim_x))
  if (!is.null(x$x_names)) {
    cat("    names           :", paste(x$x_names, collapse = ", "), "\n")
  }
  if (!is.null(x$de$embedding)) {
    cat(sprintf("  embedding (mlp)   : %d -> %d features\n",
                x$dim_x, x$de$embedding$output_dim))
  }
  cat(sprintf("  simulations       : %d\n", x$n_simulations))
  if (!is.null(x$n_dropped) && x$n_dropped > 0L) {
    total <- x$n_simulations + x$n_dropped
    cat(sprintf("    dropped         : %d of %d, non-finite output (%.1f%%)\n",
                x$n_dropped, total, 100 * x$n_dropped / total))
  }
  if (!is.null(x$de$best_val_loss) && is.finite(x$de$best_val_loss)) {
    cat(sprintf("  best val loss     : %.4f\n", x$de$best_val_loss))
  }
  if (!torch_net_alive(x$de$net)) {
    cat("  ! network unusable: a torch fit does not survive saveRDS();\n")
    cat("    save with save_npe() and reload with load_npe().\n")
  }
  cat("  -> build a posterior with posterior(fit, x_obs = ...)\n")
  invisible(x)
}
