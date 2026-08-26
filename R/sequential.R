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
#' @param simulator A function called once per parameter set, returning one
#'   simulated observation; see [nsbi_simulator].
#' @param sim_args Named list of extra arguments passed to every simulator
#'   call; see [nsbi_simulator].
#' @param x_obs The observation to target, as one numeric vector, one-row
#'   matrix or one-row data frame whose width matches the simulator's output.
#'   Sequential inference concentrates simulations around the posterior for
#'   this observation.
#' @param n_rounds Number of rounds, at least 1. Round 1 is ordinary
#'   single-round NPE.
#' @param n_simulations Simulation budget per round; either a scalar or a
#'   vector of length `n_rounds`.
#' @param density_estimator Passed to [npe()] each round.
#' @param epsilon Mass cut for the truncation: the proposal region is the
#'   `1 - epsilon` highest-probability region of the current posterior. Must
#'   be a single number strictly between 0 and 1.
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
#' prior <- prior_normal(mean = c(mu = 0, nu = 0), sd = 1)
#' simulator <- function(mu, nu) c(mu, nu) + rnorm(2, sd = 0.3)
#' fit <- npe_sequential(prior, simulator, x_obs = c(0.5, -0.5),
#'                       n_rounds = 2, n_simulations = 1000,
#'                       density_estimator = "linear_gaussian")
#' post <- posterior(fit, x_obs = c(0.5, -0.5))
#' draws <- sample(post, 1000)
#' @export
npe_sequential <- function(prior, simulator, x_obs, n_rounds = 2L,
                           n_simulations = 1000L, sim_args = list(),
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
  if (missing(x_obs)) {
    stop("`x_obs` is required: sequential NPE targets a single observation.",
         call. = FALSE)
  }
  check_x_obs(x_obs)
  # Every round simulates, so the counts are checked once, up front, rather
  # than when round r first reads one.
  n_rounds <- check_count(n_rounds, "n_rounds",
                          why = "since round 1 is the initial NPE fit")
  n_simulations <- check_counts(
    n_simulations, "n_simulations", min = 2L,
    what = "a scalar, or one budget per round")
  if (!length(n_simulations) %in% c(1L, n_rounds)) {
    stop(sprintf(paste0(
      "`n_simulations` must be length 1 or %d (`n_rounds`), not %d. ",
      "rep_len() would otherwise silently recycle a mismatched length ",
      "into per-round budgets no one asked for."),
      n_rounds, length(n_simulations)),
      call. = FALSE)
  }
  # epsilon only gets read starting in round 2 (it sets the truncation
  # threshold), but an out-of-range value should fail before round 1 spends
  # its simulation budget rather than surfacing base R's quantile() error
  # after the fact.
  epsilon <- check_prob(epsilon, "epsilon")
  n_truncation_samples <- check_count(n_truncation_samples,
                                      "n_truncation_samples")
  max_proposal_batches <- check_count(
    max_proposal_batches, "max_proposal_batches",
    why = "since one batch is one round of rejection sampling")
  if (!is.null(seed)) set.seed(seed)
  budgets <- rep_len(n_simulations, n_rounds)

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

      # keep the parameter names: they decide how the simulator is called, and
      # round 2 must call it exactly as round 1 did
      theta_new <- matrix(0, nrow = 0, ncol = prior$dim,
                          dimnames = list(NULL, prior$param_names))
      tried <- 0L
      batch <- 0L
      while (nrow(theta_new) < budgets[r] && batch < max_proposal_batches) {
        batch <- batch + 1L
        n_needed <- budgets[r] - nrow(theta_new)
        cand <- sample_prior(prior, n_needed)
        tried <- tried + nrow(cand)
        keep <- log_prob(post, cand, normalize = FALSE) >= threshold
        theta_new <- rbind(theta_new, cand[keep, , drop = FALSE])
      }
      acceptance <- nrow(theta_new) / max(tried, 1L)
      if (nrow(theta_new) == 0L) {
        # A round that adds nothing isn't "fewer simulations", it's no
        # progress at all: refitting on unchanged data would silently repeat
        # the previous round's fit and report it as round r's result.
        stop(sprintf(
          paste0("npe_sequential(): round %d accepted 0/%d proposal draws ",
                 "inside the truncated region after %d batch(es). The ",
                 "posterior estimate from round %d is too narrow relative to ",
                 "the prior for rejection sampling to find any acceptable ",
                 "proposals; try a larger `epsilon`, more ",
                 "`max_proposal_batches`, or fewer rounds."),
          r, budgets[r], batch, r - 1L), call. = FALSE)
      } else if (nrow(theta_new) < budgets[r]) {
        warning(sprintf(
          paste0("Round %d: only %d/%d proposal draws inside the truncated ",
                 "region (acceptance %.4f); continuing with fewer simulations ",
                 "this round."),
          r, nrow(theta_new), budgets[r], acceptance), call. = FALSE)
      } else {
        theta_new <- theta_new[seq_len(budgets[r]), , drop = FALSE]
      }
    }

    # `d` is unknown before round 1 has run; passing it once known keeps a
    # zero-row round (all proposals rejected by truncation) from getting a
    # 0 x 1 matrix back that can't rbind() onto x_all's real width.
    x_new <- run_simulator(simulator, theta_new, sim_args = sim_args,
                           label = sprintf("Round %d/%d", r, n_rounds),
                           d = if (r > 1L) ncol(x_all) else NULL)
    kept <- drop_failed_sims(theta_new, x_new)
    theta_new <- kept$theta
    x_new <- kept$x
    theta_all <- rbind(theta_all, theta_new)
    x_all <- rbind(x_all, x_new)
    # the simulator's output width is only known once it has run, so this is
    # the earliest point at which x_obs can be checked against it
    if (r == 1L) x_obs <- check_x_obs(x_obs, ncol(x_all))
    verbose_cat(verbose, sprintf(
      "Round %d/%d: %d new simulations (%d total), proposal acceptance %.2f\n",
      r, n_rounds, nrow(theta_new), nrow(theta_all), acceptance))

    fit <- npe(prior, theta = theta_all, x = x_all,
               density_estimator = density_estimator, verbose = verbose,
               seed = seed, ...)
    rounds[[r]] <- list(n_new = nrow(theta_new), acceptance = acceptance,
                        threshold = threshold)
  }

  fit$rounds <- rounds
  fit$x_obs <- as.numeric(x_obs)
  fit$method <- "tsnpe"
  class(fit) <- c("nsbi_snpe", class(fit))
  fit
}

#' Check the observation a sequential fit targets
#'
#' TSNPE spends its whole budget near one observation, so a mis-sized `x_obs`
#' is not a detail. Without this check `as_theta_matrix()` folds a length-3
#' vector into a 3 x 1 matrix and `resolve_x()` keeps row 1, so the run
#' truncates its proposals around a value the user never asked for, with only a
#' warning to go on. Called twice: once at the top of [npe_sequential()] for
#' shape, and once after round 1 with `dim_x`, which is the first moment the
#' simulator's output width is known.
#'
#' @param x_obs The observation passed to [npe_sequential()].
#' @param dim_x Simulator output width, or `NULL` to check shape only.
#' @return With `dim_x`, `x_obs` as a one-row matrix; otherwise `x_obs`
#'   unchanged, invisibly.
#' @keywords internal
check_x_obs <- function(x_obs, dim_x = NULL) {
  if (is.null(x_obs)) {
    stop("`x_obs` is required: sequential NPE targets a single observation.",
         call. = FALSE)
  }
  if (is.data.frame(x_obs)) x_obs <- as.matrix(x_obs)
  if (!is.numeric(x_obs) || length(x_obs) == 0L || anyNA(x_obs)) {
    stop("`x_obs` must be a non-empty numeric vector, one-row matrix or ",
         "one-row data frame, with no missing values.", call. = FALSE)
  }
  n_obs <- if (is.null(dim(x_obs))) 1L else nrow(x_obs)
  if (n_obs != 1L) {
    stop(sprintf(paste0("`x_obs` has %d rows. Sequential NPE targets one ",
                        "observation; run npe_sequential() once per ",
                        "observation, or use npe() for an amortized fit."),
                 n_obs), call. = FALSE)
  }
  if (is.null(dim_x)) return(invisible(x_obs))
  len <- if (is.null(dim(x_obs))) length(x_obs) else ncol(x_obs)
  if (len != dim_x) {
    stop(sprintf(paste0("`x_obs` has %d value(s) but the simulator returns ",
                        "%d. Sequential NPE truncates its proposals around ",
                        "`x_obs`, so a mismatch would target the wrong ",
                        "observation."), len, dim_x), call. = FALSE)
  }
  as_theta_matrix(x_obs, dim_x)
}

#' @export
print.nsbi_snpe <- function(x, ...) {
  cat("<nsbi_snpe> Sequential NPE fit (TSNPE, truncated-prior proposals)\n")
  cat(sprintf("  density estimator : %s\n", x$density_estimator))
  cat_fit_common(x, "save_npe")
  cat(sprintf("  rounds            : %d\n", length(x$rounds)))
  accs <- vapply(x$rounds, function(r) r$acceptance, numeric(1))
  cat(sprintf("  acceptance/round  : %s\n",
              paste(sprintf("%.2f", accs), collapse = ", ")))
  x_obs_display <- if (!is.null(x$x_names)) {
    paste(paste0(x$x_names, "=", signif(x$x_obs, 4)), collapse = ", ")
  } else {
    paste(signif(x$x_obs, 4), collapse = ", ")
  }
  cat(sprintf("  targeted x_obs    : %s\n", x_obs_display))
  cat("  NOT amortized: only valid at (or near) the targeted x_obs.\n")
  cat("  -> build a posterior with posterior(fit, x_obs = ...)\n")
  invisible(x)
}
