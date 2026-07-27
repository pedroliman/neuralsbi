#' MCMC over a learned likelihood
#'
#' An [nle()] fit gives an unnormalized posterior,
#' \eqn{q_\phi(x \mid \theta)\,p(\theta)}, but no way to draw from it directly.
#' `neuralsbi` samples it with a univariate slice sampler (Neal, 2003) run over
#' many chains at once.
#'
#' Slice sampling is the default in Python `sbi` for the same reason it is the
#' default here: it has no step size to tune, adapts its scale to the target as
#' it goes, and handles bounded supports without any special casing, because a
#' prior that returns `-Inf` outside its support simply shrinks the slice
#' interval.
#'
#' The vectorization runs across chains rather than across dimensions. Every
#' chain proposes its next value for the same coordinate at the same time, so
#' one step costs one batched call to the density estimator instead of
#' `n_chains` separate forward passes. With a neural likelihood that difference
#' is the whole running time.
#'
#' @references Neal, R. M. (2003). Slice sampling. *The Annals of Statistics*
#'   31(3), 705-767.
#'
#' @name nsbi_mcmc
#' @keywords internal
NULL

#' Vectorized univariate slice sampler
#'
#' @param log_prob_fn Vectorized log-density: takes an `m x dim` matrix, returns
#'   `m` log-densities (`-Inf` allowed).
#' @param init `n_chains x dim` matrix of starting points, all with finite
#'   log-density.
#' @param n_draws Number of retained draws in total, across all chains.
#' @param warmup Steps discarded at the start of each chain.
#' @param thin Keep one draw in `thin`.
#' @param width Initial slice width per dimension (recycled). Adapted during
#'   warmup and then held fixed.
#' @param max_steps Cap on stepping-out expansions per coordinate.
#' @param verbose Report progress.
#' @return A list with `draws` (`n_draws x dim`), `chains` (`n_kept x n_chains x
#'   dim` array) and `n_evals`.
#' @keywords internal
slice_sample <- function(log_prob_fn, init, n_draws, warmup = 200L,
                         thin = 10L, width = 1, max_steps = 100L,
                         verbose = FALSE) {
  with_nsbi_progress(slice_sample_run(log_prob_fn, init, n_draws, warmup, thin,
                                      width, max_steps, verbose))
}

#' The loop behind [slice_sample()], split out so the progress context wraps it
#' @keywords internal
slice_sample_run <- function(log_prob_fn, init, n_draws, warmup, thin, width,
                             max_steps, verbose) {
  init <- as_theta_matrix(init)
  n_chains <- nrow(init)
  dim <- ncol(init)
  width <- rep_len(as.numeric(width), dim)

  n_kept <- ceiling(n_draws / n_chains)
  n_iter <- warmup + n_kept * thin

  state <- init
  state_lp <- log_prob_fn(state)
  if (any(!is.finite(state_lp))) {
    stop("Some MCMC starting points have zero posterior density. ",
         "This is an initialization failure, not a sampling one.", call. = FALSE)
  }
  n_evals <- n_chains

  kept <- array(0, dim = c(n_kept, n_chains, dim))
  keep_i <- 0L
  chain_ids <- seq_len(n_chains)

  p <- nsbi_progressor(steps = n_iter, label = "Sampling")

  for (iter in seq_len(n_iter)) {
    for (d in seq_len(dim)) {
      x0 <- state[, d]
      # The slice level: log(u) below the current density, per chain.
      level <- state_lp - stats::rexp(n_chains)

      # -- stepping out ---------------------------------------------------
      # Random placement of the initial interval keeps the sampler reversible.
      lo <- x0 - width[d] * stats::runif(n_chains)
      hi <- lo + width[d]

      lo_open <- rep(TRUE, n_chains)
      for (s in seq_len(max_steps)) {
        active <- which(lo_open)
        if (!length(active)) break
        cand <- state[active, , drop = FALSE]
        cand[, d] <- lo[active]
        lp <- log_prob_fn(cand)
        n_evals <- n_evals + length(active)
        grow <- lp > level[active]
        lo[active[grow]] <- lo[active[grow]] - width[d]
        lo_open[active[!grow]] <- FALSE
      }
      hi_open <- rep(TRUE, n_chains)
      for (s in seq_len(max_steps)) {
        active <- which(hi_open)
        if (!length(active)) break
        cand <- state[active, , drop = FALSE]
        cand[, d] <- hi[active]
        lp <- log_prob_fn(cand)
        n_evals <- n_evals + length(active)
        grow <- lp > level[active]
        hi[active[grow]] <- hi[active[grow]] + width[d]
        hi_open[active[!grow]] <- FALSE
      }

      # -- shrinkage ------------------------------------------------------
      pending <- chain_ids
      while (length(pending)) {
        prop <- lo[pending] + stats::runif(length(pending)) * (hi[pending] - lo[pending])
        cand <- state[pending, , drop = FALSE]
        cand[, d] <- prop
        lp <- log_prob_fn(cand)
        n_evals <- n_evals + length(pending)

        accept <- lp > level[pending]
        acc_id <- pending[accept]
        if (length(acc_id)) {
          state[acc_id, d] <- prop[accept]
          state_lp[acc_id] <- lp[accept]
        }
        # A rejected proposal becomes the new interval edge on its own side.
        rej <- pending[!accept]
        if (length(rej)) {
          pr <- prop[!accept]
          left <- pr < x0[rej]
          lo[rej[left]] <- pr[left]
          hi[rej[!left]] <- pr[!left]
          # Guard against an interval collapsing to a point on a flat or
          # numerically degenerate target. Such a chain stops shrinking and
          # keeps its current value rather than looping forever.
          stuck <- rej[(hi[rej] - lo[rej]) <= .Machine$double.eps * 8]
          if (length(stuck)) rej <- setdiff(rej, stuck)
        }
        pending <- rej
      }
    }

    if (iter <= warmup) {
      # Adapt the slice width to the target's own scale. A width set from the
      # prior is badly wrong once the posterior is much narrower -- which is
      # exactly what many independent observations produce -- and every unit of
      # mismatch is paid for in stepping-out or shrinkage evaluations. The
      # across-chain spread is a free estimate of that scale. Adaptation stops
      # at the end of warmup, so the retained draws come from a fixed kernel.
      spread <- apply(state, 2, stats::sd)
      ok <- is.finite(spread) & spread > 0
      width[ok] <- 0.5 * width[ok] + 0.5 * (2 * spread[ok])
    } else if ((iter - warmup) %% thin == 0L) {
      keep_i <- keep_i + 1L
      kept[keep_i, , ] <- state
    }
    p(1, done = iter == n_iter)
  }

  # Interleave chains so a truncated result still spans all of them.
  draws <- matrix(aperm(kept, c(2, 1, 3)), ncol = dim)
  draws <- draws[seq_len(min(n_draws, nrow(draws))), , drop = FALSE]
  list(draws = draws, chains = kept, n_evals = n_evals)
}

#' Starting points for the chains
#'
#' `"resample"` (the default, and `sbi`'s) draws a large pool from the prior,
#' weights it by the posterior density and resamples without replacement: a
#' sampling-importance-resampling start that puts the chains where the mass is,
#' which matters because slice sampling has no adaptation phase to rescue a bad
#' start. `"proposal"` just takes prior draws.
#'
#' @keywords internal
mcmc_init <- function(prior, log_prob_fn, n_chains,
                      strategy = c("resample", "proposal"),
                      n_pool = 1000L) {
  strategy <- match.arg(strategy)
  if (strategy == "proposal") {
    for (attempt in seq_len(20L)) {
      cand <- sample_prior(prior, n_chains)
      if (all(is.finite(log_prob_fn(cand)))) return(cand)
    }
    stop("Could not find prior draws with finite posterior density after 20 ",
         "attempts. The surrogate likelihood may be degenerate.", call. = FALSE)
  }

  pool <- sample_prior(prior, max(n_pool, n_chains))
  lp <- log_prob_fn(pool)
  ok <- which(is.finite(lp))
  if (!length(ok)) {
    stop("No prior draw has finite posterior density, so the chains cannot be ",
         "started. Check the prior and the fitted likelihood.", call. = FALSE)
  }
  if (length(ok) <= n_chains) {
    return(pool[rep_len(ok, n_chains), , drop = FALSE])
  }
  w <- exp(lp[ok] - max(lp[ok]))
  idx <- ok[sample.int(length(ok), n_chains, replace = FALSE, prob = w)]
  pool[idx, , drop = FALSE]
}

#' Split-Rhat and bulk effective sample size
#'
#' The standard rank-free versions from Vehtari et al. (2021), computed on the
#' split chains. Implemented here rather than taken from \pkg{posterior} to keep
#' the dependency surface where it is; the test suite cross-checks against
#' \pkg{posterior} when that package happens to be installed.
#'
#' @param chains A `n_iter x n_chains x dim` array.
#' @return A data frame with one row per parameter: `rhat` and `ess_bulk`.
#' @references Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and
#'   Burkner, P.-C. (2021). Rank-normalization, folding, and localization.
#'   *Bayesian Analysis* 16(2), 667-718.
#' @keywords internal
mcmc_diagnostics <- function(chains) {
  n_iter <- dim(chains)[1]
  n_chains <- dim(chains)[2]
  dimension <- dim(chains)[3]
  out <- data.frame(rhat = rep(NA_real_, dimension),
                    ess_bulk = rep(NA_real_, dimension))
  if (n_iter < 4L) return(out)

  half <- floor(n_iter / 2)
  for (j in seq_len(dimension)) {
    # Split each chain in two, so a chain that drifts is caught.
    m <- cbind(chains[seq_len(half), , j, drop = FALSE][, , 1],
               chains[n_iter - half + seq_len(half), , j, drop = FALSE][, , 1])
    m <- matrix(m, nrow = half)
    out$rhat[j] <- split_rhat(m)
    out$ess_bulk[j] <- bulk_ess(m)
  }
  out
}

#' @keywords internal
split_rhat <- function(m) {
  n <- nrow(m)
  k <- ncol(m)
  if (n < 2L || k < 2L) return(NA_real_)
  chain_means <- colMeans(m)
  chain_vars <- apply(m, 2, stats::var)
  if (any(!is.finite(chain_vars))) return(NA_real_)
  W <- mean(chain_vars)
  B <- n * stats::var(chain_means)
  if (W <= 0) return(NA_real_)
  var_hat <- ((n - 1) / n) * W + B / n
  sqrt(var_hat / W)
}

#' @keywords internal
bulk_ess <- function(m) {
  n <- nrow(m)
  k <- ncol(m)
  if (n < 4L) return(NA_real_)
  # Rank-normalize, which is what makes this robust to heavy tails.
  r <- matrix(rank(as.vector(m), ties.method = "average"), nrow = n)
  z <- stats::qnorm((r - 0.375) / (n * k + 0.25))

  chain_vars <- apply(z, 2, stats::var)
  W <- mean(chain_vars)
  if (!is.finite(W) || W <= 0) return(NA_real_)
  B <- n * stats::var(colMeans(z))
  var_hat <- ((n - 1) / n) * W + B / n

  # Average autocorrelation across chains, via the FFT. autocov() divides by n,
  # so the n/(n-1) factor puts it on the same footing as the unbiased chain
  # variances in W.
  acov <- rowMeans(vapply(seq_len(k), function(c) autocov(z[, c], n), numeric(n)))
  rho <- 1 - (W - acov * (n / (n - 1))) / var_hat

  # Geyer's initial positive sequence: sum paired autocorrelations while the
  # pair sum stays positive, then enforce the monotone decrease the theory
  # guarantees but a finite sample does not.
  pairs <- numeric(0)
  t <- 0L                     # the first pair is lag 0 + lag 1
  while (t + 2L <= n) {
    pair <- rho[t + 1L] + rho[t + 2L]
    if (!is.finite(pair) || pair <= 0) break
    pairs <- c(pairs, pair)
    t <- t + 2L
  }
  if (length(pairs) > 1L) pairs <- cummin(pairs)
  tau <- -1 + 2 * sum(pairs)
  if (!is.finite(tau) || tau < 1) tau <- 1
  min(n * k / tau, n * k * log10(n * k))
}

#' Autocovariance of a single chain at lags 0..n-1, via the FFT
#' @keywords internal
autocov <- function(y, n) {
  y <- y - mean(y)
  nfft <- stats::nextn(2 * n, 2)
  f <- stats::fft(c(y, rep(0, nfft - n)))
  ac <- Re(stats::fft(f * Conj(f), inverse = TRUE)) / nfft
  ac[seq_len(n)] / n
}
