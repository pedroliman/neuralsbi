#' Sequential NPE with truncated-prior proposals (TSNPE)
#'
#' Multi-round NPE targeting a single observation `x_obs`. Single-round [npe()]
#' spends its simulation budget across the whole prior; when only one
#' observation matters, most of those simulations land in regions the posterior
#' never visits. `npe_sequential()` implements truncated sequential NPE (TSNPE,
#' Deistler et al. 2022): after each round the prior is truncated to the
#' highest-probability region of the current posterior estimate, and the next
#' round's parameters are drawn from that truncated prior. Because every
#' proposal is proportional to the prior on its support, the standard NPE loss
#' stays valid -- no importance-weight or atomic correction is needed, which is
#' what makes TSNPE the simplest correct sequential scheme.
#'
#' The rounds accumulate: each round's estimator is trained on all simulations
#' so far. The final fit is returned as an `nsbi_npe` (subclass `nsbi_snpe`)
#' and works with [posterior()], [sample()] and the diagnostics, but unlike
#' single-round NPE it is *not* amortized: it is only trustworthy at (or very
#' near) `x_obs`.
#'
#' Proposal draws are obtained by rejection: prior candidates are kept when
#' their posterior log-density clears the `epsilon`-quantile threshold of the
#' current posterior's own draws. If the posterior is much narrower than the
#' prior the acceptance rate falls; the round then stops after
#' `max_proposal_batches` batches and continues with the draws it has,
#' with a warning.
#'
#' @param prior An `nsbi_prior` (see [prior_uniform()], [prior_normal()]).
#' @param simulator A function mapping an `n x dim` matrix of parameters to an
#'   `n x d` matrix of simulated data.
#' @param x_obs The observation to target. Sequential inference concentrates
#'   simulations around the posterior for this observation.
#' @param n_rounds Number of rounds. Round 1 is ordinary single-round NPE.
#' @param n_simulations Simulation budget per round; either a scalar or a
#'   vector of length `n_rounds`.
#' @param density_estimator Passed to [npe()] each round.
#' @param epsilon Mass cut for the truncation: the proposal region is the
#'   `1 - epsilon` highest-probability region of the current posterior.
#' @param n_truncation_samples Posterior draws used to locate the truncation
#'   threshold each round.
#' @param max_proposal_batches Cap on rejection-sampling batches per round.
#' @param seed Optional integer seed for reproducibility.
#' @param verbose Print per-round progress.
#' @param ... Passed to [npe()] (estimator and training settings).
#'
#' @return An object of class `c("nsbi_snpe", "nsbi_npe")` with a `rounds`
#'   field recording per-round budgets, acceptance rates, and thresholds.
#'
#' @references Deistler, Goncalves & Macke (2022), "Truncated proposals for
#'   scalable and hassle-free simulation-based inference", NeurIPS.
#'   \doi{10.48550/arXiv.2210.04815}
#'
#' @examples
#' prior <- prior_normal(mean = c(0, 0), sd = 1)
#' simulator <- function(theta) theta + matrix(rnorm(length(theta), sd = 0.3),
#'                                             nrow = nrow(theta))
#' fit <- npe_sequential(prior, simulator, x_obs = c(0.5, -0.5),
#'                       n_rounds = 2, n_simulations = 1000,
#'                       density_estimator = "linear_gaussian")
#' post <- posterior(fit, x_obs = c(0.5, -0.5))
#' draws <- sample(post, 1000)
#' @export
npe_sequential <- function(prior, simulator, x_obs, n_rounds = 2L,
                           n_simulations = 1000L,
                           density_estimator = c("maf", "mdn", "nsf",
                                                 "linear_gaussian"),
                           epsilon = 1e-4, n_truncation_samples = 5000L,
                           max_proposal_batches = 200L,
                           seed = NULL, verbose = FALSE, ...) {
  stopifnot(inherits(prior, "nsbi_prior"))
  if (!is.function(simulator)) {
    stop("`simulator` must be a function; sequential NPE has to simulate ",
         "each round.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  n_rounds <- as.integer(n_rounds)
  budgets <- rep_len(as.integer(n_simulations), n_rounds)

  theta_all <- NULL
  x_all <- NULL
  fit <- NULL
  rounds <- vector("list", n_rounds)

  for (r in seq_len(n_rounds)) {
    if (r == 1L) {
      theta_new <- sample_prior(prior, budgets[r])
      acceptance <- 1
      threshold <- -Inf
    } else {
      post <- posterior(fit, x_obs = x_obs)
      ref <- sample.nsbi_posterior(post, n = n_truncation_samples)
      lp_ref <- log_prob(post, ref, normalize = FALSE)
      threshold <- stats::quantile(lp_ref, probs = epsilon, names = FALSE)

      theta_new <- matrix(0, nrow = 0, ncol = prior$dim)
      tried <- 0L
      batch <- 0L
      while (nrow(theta_new) < budgets[r] && batch < max_proposal_batches) {
        batch <- batch + 1L
        cand <- sample_prior(prior, budgets[r])
        tried <- tried + nrow(cand)
        keep <- log_prob(post, cand, normalize = FALSE) >= threshold
        theta_new <- rbind(theta_new, cand[keep, , drop = FALSE])
      }
      acceptance <- nrow(theta_new) / max(tried, 1L)
      if (nrow(theta_new) < budgets[r]) {
        warning(sprintf(
          paste0("Round %d: only %d/%d proposal draws inside the truncated ",
                 "region (acceptance %.4f); continuing with fewer simulations ",
                 "this round."),
          r, nrow(theta_new), budgets[r], acceptance), call. = FALSE)
      } else {
        theta_new <- theta_new[seq_len(budgets[r]), , drop = FALSE]
      }
    }

    x_new <- as_theta_matrix(simulator(theta_new))
    if (nrow(x_new) != nrow(theta_new)) {
      stop("Simulator must return one row of output per row of theta.",
           call. = FALSE)
    }
    theta_all <- rbind(theta_all, theta_new)
    x_all <- rbind(x_all, x_new)
    verbose_cat(verbose, sprintf(
      "Round %d/%d: %d new simulations (%d total), proposal acceptance %.2f\n",
      r, n_rounds, nrow(theta_new), nrow(theta_all), acceptance))

    fit <- npe(prior, theta = theta_all, x = x_all,
               density_estimator = density_estimator, verbose = verbose, ...)
    rounds[[r]] <- list(n_new = nrow(theta_new), acceptance = acceptance,
                        threshold = threshold)
  }

  fit$rounds <- rounds
  fit$x_obs <- as.numeric(as_theta_matrix(x_obs, fit$dim_x))
  fit$method <- "tsnpe"
  class(fit) <- c("nsbi_snpe", class(fit))
  fit
}

#' @export
print.nsbi_snpe <- function(x, ...) {
  cat("<nsbi_snpe> Sequential NPE fit (TSNPE, truncated-prior proposals)\n")
  cat(sprintf("  density estimator : %s\n", x$density_estimator))
  cat(sprintf("  rounds            : %d\n", length(x$rounds)))
  cat(sprintf("  simulations       : %d\n", x$n_simulations))
  accs <- vapply(x$rounds, function(r) r$acceptance, numeric(1))
  cat(sprintf("  acceptance/round  : %s\n",
              paste(sprintf("%.2f", accs), collapse = ", ")))
  cat(sprintf("  targeted x_obs    : %s\n",
              paste(signif(x$x_obs, 4), collapse = ", ")))
  cat("  NOT amortized: only valid at (or near) the targeted x_obs.\n")
  cat("  -> build a posterior with posterior(fit, x_obs = ...)\n")
  invisible(x)
}
