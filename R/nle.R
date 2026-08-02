#' Neural Likelihood Estimation (NLE)
#'
#' `nle()` trains a conditional density estimator on the *other* factorization
#' of the joint: instead of the posterior \eqn{p(\theta \mid x)} that [npe()]
#' learns, it learns a surrogate likelihood \eqn{q_\phi(x \mid \theta)}. The
#' posterior then follows from Bayes' rule,
#' \eqn{p(\theta \mid x) \propto q_\phi(x \mid \theta)\,p(\theta)}, and is
#' sampled with MCMC by [posterior()].
#'
#' @section When NLE beats NPE:
#'
#' An NPE fit learns the posterior for one fixed data dimension, chosen at
#' training time. If the observation is \eqn{n} independent trials from the same
#' parameter, NPE must be retrained for every \eqn{n}, or handed summary
#' statistics that throw information away. NLE learns the density of a *single*
#' trial, so the log-likelihood of \eqn{n} trials is a sum of \eqn{n}
#' evaluations:
#'
#' \deqn{\log p(x_1, \ldots, x_n \mid \theta) = \sum_{i=1}^{n} \log q_\phi(x_i \mid \theta).}
#'
#' Train once, then condition on 50 trials or 5000 without touching the network.
#' The learned likelihood is also a plain differentiable function of
#' \eqn{\theta}, so it can be embedded in a larger model written by hand --
#' see [stan_code()].
#'
#' The trade-off is real and worth stating: posterior draws now cost an MCMC
#' run rather than a forward pass, and for a single fixed observation with
#' high-dimensional data NPE is usually the better choice.
#'
#' @inheritParams npe
#' @param density_estimator One of `"maf"` (Masked Autoregressive Flow, needs
#'   `torch`; the default), `"mdn"` (Mixture Density Network, needs `torch`),
#'   `"nsf"` (Neural Spline Flow, needs `torch`), or `"linear_gaussian"`
#'   (closed-form baseline, no `torch`), or a function `function(theta, x)`
#'   returning a fitted estimator. Note the estimator sees the roles swapped:
#'   its target is `x` and it conditions on `theta`.
#' @param n_transforms MAF/NSF setting: number of stacked autoregressive
#'   transforms.
#' @param n_components,hidden MDN settings: number of mixture components and a
#'   vector of hidden-layer widths.
#' @param n_bins,tail_bound NSF settings.
#'
#' @return An object of class `nsbi_nle`. Evaluate the surrogate likelihood with
#'   [log_lik()], turn it into a posterior with [posterior()], or export it to
#'   Stan with [stan_code()].
#'
#' @seealso [log_lik()] and [likelihood_fn()] to evaluate the surrogate,
#'   [posterior()] to sample it, [stan_code()] to hand it to Stan.
#'
#' @references Papamakarios, G., Sterratt, D. and Murray, I. (2019). Sequential
#'   Neural Likelihood: Fast Likelihood-free Inference with Autoregressive
#'   Flows. *AISTATS*.
#'
#' @examples
#' # One noisy measurement per simulator call; the observation is 200 of them.
#' prior <- prior_uniform(c(mu = -3, log_sigma = -1),
#'                        c(mu =  3, log_sigma =  1))
#' simulator <- function(mu, log_sigma) c(y = rnorm(1, mu, exp(log_sigma)))
#'
#' fit <- nle(prior, simulator, n_simulations = 2000,
#'            density_estimator = "linear_gaussian")
#'
#' x_obs <- matrix(rnorm(200, mean = 1, sd = 0.5), ncol = 1)
#' log_lik(fit, theta = c(1, log(0.5)), x = x_obs)
#' @export
nle <- function(prior, simulator = NULL, n_simulations = 1000,
                sim_args = list(), theta = NULL, x = NULL,
                density_estimator = c("maf", "mdn", "nsf", "linear_gaussian"),
                n_components = 10L, n_transforms = 5L, hidden = c(50L, 50L),
                n_bins = 10L, tail_bound = 3,
                max_epochs = 2000L, batch_size = 200L, lr = 5e-4,
                validation_fraction = 0.1, patience = 20L,
                n_restarts = 1L, clip_grad_norm = 5,
                standardize = TRUE, seed = NULL, verbose = FALSE) {
  # See npe(): everything here is checked before the simulator runs.
  if (!is.function(density_estimator)) {
    density_estimator <- match.arg(density_estimator)
  }
  check_prior(prior)
  check_architecture(n_components, n_transforms, hidden, n_bins, tail_bound)
  check_train_controls(max_epochs, batch_size, lr, validation_fraction,
                       patience, n_restarts, clip_grad_norm)

  prep <- prepare_simulations(prior, simulator, n_simulations, sim_args,
                              theta, x, standardize, seed, verbose)

  # The estimator contract is q(target | condition) and makes no assumption
  # about which of the two is a parameter. NLE is NPE with the roles swapped:
  # the target is the data, the conditioning variable is the parameter.
  de <- fit_density_estimator(
    density_estimator, theta_z = prep$x_z, x_z = prep$theta_z,
    n_components = n_components, n_transforms = n_transforms,
    hidden = hidden, n_bins = n_bins, tail_bound = tail_bound,
    embedding_net = NULL, max_epochs = max_epochs,
    batch_size = batch_size, lr = lr, validation_fraction = validation_fraction,
    patience = patience, n_restarts = n_restarts,
    clip_grad_norm = clip_grad_norm, seed = seed, verbose = verbose
  )

  structure(
    list(
      de = de,
      prior = prior,
      std_theta = prep$std_theta,
      std_x = prep$std_x,
      dim_theta = prior$dim,
      dim_x = ncol(prep$x),
      param_names = prep$param_names,
      x_names = prep$x_names,
      n_simulations = nrow(prep$theta),
      n_dropped = prep$n_dropped,
      density_estimator = if (is.character(density_estimator))
        density_estimator[1] else "custom"
    ),
    class = "nsbi_nle"
  )
}

#' @export
print.nsbi_nle <- function(x, ...) {
  cat("<nsbi_nle> Neural Likelihood Estimation fit\n")
  cat(sprintf("  density estimator : %s  (learns q(x | theta))\n",
              x$density_estimator))
  cat(sprintf("  parameters (dim)  : %d\n", x$dim_theta))
  if (!is.null(x$param_names)) {
    cat("    names           :", paste(x$param_names, collapse = ", "), "\n")
  }
  cat(sprintf("  data (dim)        : %d  per observation\n", x$dim_x))
  if (!is.null(x$x_names)) {
    cat("    names           :", paste(x$x_names, collapse = ", "), "\n")
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
    cat("    save with save_nle() and reload with load_nle().\n")
  }
  cat("  -> log_lik(fit, theta, x), posterior(fit, x_obs = ...), stan_code(fit)\n")
  invisible(x)
}
