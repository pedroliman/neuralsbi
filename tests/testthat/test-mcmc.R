test_that("the slice sampler recovers a standard normal", {
  set.seed(1)
  lp <- function(theta) rowSums(stats::dnorm(theta, log = TRUE))
  init <- matrix(stats::rnorm(8 * 2), ncol = 2)

  res <- slice_sample(lp, init, n_draws = 4000, warmup = 200, thin = 5)

  expect_equal(ncol(res$draws), 2L)
  expect_equal(nrow(res$draws), 4000L)
  expect_equal(colMeans(res$draws), c(0, 0), tolerance = 0.1)
  expect_equal(apply(res$draws, 2, stats::sd), c(1, 1), tolerance = 0.1)
})

test_that("the slice sampler recovers a correlated Gaussian", {
  set.seed(2)
  Sigma <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
  P <- solve(Sigma)
  lp <- function(theta) -0.5 * rowSums((theta %*% P) * theta)
  init <- matrix(stats::rnorm(10 * 2, sd = 0.5), ncol = 2)

  res <- slice_sample(lp, init, n_draws = 6000, warmup = 300, thin = 5)

  expect_equal(colMeans(res$draws), c(0, 0), tolerance = 0.12)
  expect_equal(stats::cor(res$draws)[1, 2], 0.8, tolerance = 0.06)
})

test_that("the sampler never leaves a bounded support", {
  set.seed(3)
  lp <- function(theta) {
    ifelse(theta[, 1] >= -1 & theta[, 1] <= 1, 0, -Inf)
  }
  init <- matrix(stats::runif(6, -0.5, 0.5), ncol = 1)

  res <- slice_sample(lp, init, n_draws = 3000, warmup = 100, thin = 2,
                      width = 0.5)

  expect_true(all(res$draws >= -1 & res$draws <= 1))
  # A uniform on [-1, 1] has mean 0 and variance 1/3.
  expect_equal(mean(res$draws), 0, tolerance = 0.06)
  expect_equal(stats::var(as.numeric(res$draws)), 1 / 3, tolerance = 0.05)
})

test_that("one sweep of the kernel leaves the target alone", {
  # The strongest statement available about a sampler: start the chains at
  # exact draws from the target, take a single sweep, and the result must still
  # be exact draws from the target. Recovering the moments of a normal, as the
  # tests above do, is a weaker check -- it passes for a kernel that is slightly
  # wrong in a direction the mean and sd do not see.
  #
  # This matters here because the sweep is not the textbook one. Stepping out
  # carries several of an edge's future positions in one call and shrinkage
  # several of a chain's future proposals, on the argument that both sequences
  # are determined before any density is known. If that argument is wrong
  # anywhere, invariance is what breaks.
  #
  # One chain per draw, one sweep, no warmup, so the draws are independent and
  # the KS p-value is calibrated. No warmup also means no width adaptation, so
  # this is a fixed kernel rather than an adapting one.
  set.seed(12)
  n <- 600L

  normal_lp <- function(theta) stats::dnorm(theta[, 1], log = TRUE)
  from_normal <- slice_sample(normal_lp, matrix(stats::rnorm(n), ncol = 1),
                              n_draws = n, warmup = 0L, thin = 1L, width = 1)
  expect_gt(stats::ks.test(as.numeric(from_normal$draws), "pnorm")$p.value, 0.01)

  # A skewed, bounded-below target: a kernel that treats the two sides of the
  # current point differently shows up here and not on a symmetric one.
  gamma_lp <- function(theta) stats::dgamma(theta[, 1], shape = 2, log = TRUE)
  from_gamma <- slice_sample(gamma_lp, matrix(stats::rgamma(n, 2), ncol = 1),
                             n_draws = n, warmup = 0L, thin = 1L, width = 0.3)
  expect_gt(stats::ks.test(as.numeric(from_gamma$draws), "pgamma",
                           shape = 2)$p.value, 0.01)

  # A width far narrower than the target forces the stepping-out loop deep,
  # which is the path that batches several expansions into one call.
  wide_lp <- function(theta) stats::dnorm(theta[, 1], sd = 10, log = TRUE)
  from_wide <- slice_sample(wide_lp, matrix(stats::rnorm(n, sd = 10), ncol = 1),
                            n_draws = n, warmup = 0L, thin = 1L, width = 0.2)
  expect_gt(stats::ks.test(as.numeric(from_wide$draws), "pnorm",
                           sd = 10)$p.value, 0.01)
})

test_that("warmup and thinning give the requested number of draws", {
  set.seed(4)
  lp <- function(theta) rowSums(stats::dnorm(theta, log = TRUE))
  init <- matrix(0, nrow = 4, ncol = 1)

  res <- slice_sample(lp, init, n_draws = 100, warmup = 10, thin = 3)

  expect_equal(nrow(res$draws), 100L)
  expect_equal(dim(res$chains), c(25L, 4L, 1L))
})

test_that("a starting point with zero density is refused", {
  lp <- function(theta) rep(-Inf, nrow(theta))
  expect_error(
    slice_sample(lp, matrix(0, 2, 1), n_draws = 10, warmup = 1, thin = 1),
    "zero posterior density"
  )
})

test_that("a non-finite or non-positive width is refused, not left to crash the stepping-out loop", {
  # GitHub #164: width used to reach `stats::runif()` and the array indexing
  # inside the stepping-out loop unchecked, so `width = NA` failed as "NAs are
  # not allowed in subscripted assignments" -- several frames deep, naming
  # neither `width` nor the sampler.
  lp <- function(theta) rowSums(stats::dnorm(theta, log = TRUE))
  init <- matrix(0, nrow = 4, ncol = 1)

  expect_error(
    slice_sample(lp, init, n_draws = 10, warmup = 1, thin = 1, width = NA),
    "`width`"
  )
  expect_error(
    slice_sample(lp, init, n_draws = 10, warmup = 1, thin = 1, width = NaN),
    "`width`"
  )
  expect_error(
    slice_sample(lp, init, n_draws = 10, warmup = 1, thin = 1, width = Inf),
    "`width`"
  )
  expect_error(
    slice_sample(lp, init, n_draws = 10, warmup = 1, thin = 1, width = 0),
    "`width`"
  )
  expect_error(
    slice_sample(lp, init, n_draws = 10, warmup = 1, thin = 1, width = -1),
    "`width`"
  )

  # dim = 2: one bad entry among several should still be caught after
  # rep_len() recycling.
  init2 <- matrix(0, nrow = 4, ncol = 2)
  lp2 <- function(theta) rowSums(stats::dnorm(theta, log = TRUE))
  expect_error(
    slice_sample(lp2, init2, n_draws = 10, warmup = 1, thin = 1,
                 width = c(1, NA)),
    "`width`"
  )
})

test_that("a non-finite or non-positive max_steps is refused, not left to loop on a bad count", {
  # GitHub #190: the NEWS entry for #164/#167 flagged max_steps and n_pool as
  # the same gap width had before that fix -- steps_left <- max_steps drives
  # the stepping-out loop's `while (steps_left > 0L)` with no validation.
  lp <- function(theta) rowSums(stats::dnorm(theta, log = TRUE))
  init <- matrix(0, nrow = 4, ncol = 1)

  expect_error(
    slice_sample(lp, init, n_draws = 10, warmup = 1, thin = 1, max_steps = NA),
    "`max_steps`"
  )
  expect_error(
    slice_sample(lp, init, n_draws = 10, warmup = 1, thin = 1, max_steps = 0),
    "`max_steps`"
  )
  expect_error(
    slice_sample(lp, init, n_draws = 10, warmup = 1, thin = 1, max_steps = -1),
    "`max_steps`"
  )
  expect_error(
    slice_sample(lp, init, n_draws = 10, warmup = 1, thin = 1, max_steps = 1.5),
    "`max_steps`"
  )
})

test_that("mcmc_init puts chains where the mass is", {
  set.seed(5)
  prior <- prior_uniform(low = -10, high = 10)
  # Density concentrated near 5: the resampling start should find it, the
  # plain prior start should not.
  lp <- function(theta) stats::dnorm(theta[, 1], mean = 5, sd = 0.5, log = TRUE)

  resampled <- mcmc_init(prior, lp, n_chains = 20, strategy = "resample")
  proposed <- mcmc_init(prior, lp, n_chains = 20, strategy = "proposal")

  expect_equal(dim(resampled), c(20L, 1L))
  expect_lt(mean(abs(resampled - 5)), mean(abs(proposed - 5)))
})

test_that("resampling survives a likelihood peaked enough to underflow", {
  # The case this broke on: a few thousand independent observations make the
  # log-density spread across prior draws run to thousands, so exponentiating
  # the weights before sampling leaves almost all of them at exactly zero.
  set.seed(11)
  prior <- prior_uniform(low = -10, high = 10)
  lp <- function(theta) 5000 * stats::dnorm(theta[, 1], 5, 0.5, log = TRUE)

  init <- mcmc_init(prior, lp, n_chains = 20, strategy = "resample")

  expect_equal(dim(init), c(20L, 1L))
  expect_true(all(is.finite(lp(init))))
  expect_lt(mean(abs(init - 5)), 2)
})

test_that("mcmc_init reports an impossible start rather than looping", {
  prior <- prior_uniform(low = 0, high = 1)
  expect_error(
    mcmc_init(prior, function(theta) rep(-Inf, nrow(theta)), n_chains = 2),
    "finite posterior density"
  )
})

test_that("mcmc_init proposal strategy gives up after 20 attempts", {
  prior <- prior_uniform(low = 0, high = 1)
  expect_error(
    mcmc_init(prior, function(theta) rep(-Inf, nrow(theta)), n_chains = 2,
              strategy = "proposal"),
    "after 20 attempts"
  )
})

test_that("mcmc_init proposal strategy pools across attempts rather than requiring one lucky batch", {
  # Finite on 65% of the prior's support: entirely ordinary, not degenerate.
  # A single batch of n_chains = 20 draws is jointly finite with probability
  # 0.65^20 ~ 3e-4, so requiring one batch to succeed whole used to fail here
  # almost every time across 20 retries.
  set.seed(1)
  prior <- prior_uniform(low = c(-10, -10), high = c(10, 10))
  lp <- function(theta) ifelse(theta[, 1] > 3, -Inf, 0)

  init <- mcmc_init(prior, lp, n_chains = 20, strategy = "proposal")

  expect_equal(dim(init), c(20L, 2L))
  expect_true(all(is.finite(lp(init))))
})

test_that("mcmc_init proposal strategy reports a low acceptance rate, not degeneracy", {
  # Acceptance is about 1%, so the 20-attempt budget (n_pool <= n_chains keeps
  # each batch at n_chains draws) collects roughly 200 finite draws in
  # expectation -- far short of the 1000 needed, but nowhere near zero, so
  # this is a low-acceptance failure, not a degenerate one.
  set.seed(2)
  prior <- prior_uniform(low = 0, high = 1)
  lp <- function(theta) ifelse(theta[, 1] < 0.01, 0, -Inf)

  err <- tryCatch(
    mcmc_init(prior, lp, n_chains = 1000, strategy = "proposal", n_pool = 1),
    error = function(e) e
  )
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "low acceptance rate")
  expect_no_match(conditionMessage(err), "may be degenerate")
})

test_that("mcmc_init resample strategy pools across attempts instead of recycling", {
  # GitHub #182: a single pool of n_pool draws can land fewer than n_chains of
  # them in the posterior's support. The old code padded the shortfall with
  # rep_len(ok, n_chains), which starts several chains from identical points
  # -- exactly what split-Rhat and bulk ESS assume never happens. Finite only
  # for theta[, 1] > 9.5 out of a [-10, 10] prior is a narrow sliver: one pool
  # of 1000 draws lands only a handful there, far fewer than n_chains = 20.
  set.seed(1)
  prior <- prior_uniform(low = c(-10, -10), high = c(10, 10))
  lp <- function(theta) ifelse(theta[, 1] > 9.5, 0, -Inf)

  init <- mcmc_init(prior, lp, n_chains = 20, strategy = "resample", n_pool = 1000)

  expect_equal(dim(init), c(20L, 2L))
  expect_true(all(is.finite(lp(init))))
  expect_equal(nrow(unique(init)), 20L)
})

test_that("mcmc_init resample strategy reports a low acceptance rate, not degeneracy", {
  # Acceptance is about 1%, so a single pool of n_pool = n_chains draws rarely
  # finds n_chains finite ones, but the 20-attempt budget collects roughly 200
  # in expectation -- far short of the 1000 needed, but nowhere near zero.
  set.seed(2)
  prior <- prior_uniform(low = 0, high = 1)
  lp <- function(theta) ifelse(theta[, 1] < 0.01, 0, -Inf)

  err <- tryCatch(
    mcmc_init(prior, lp, n_chains = 1000, strategy = "resample", n_pool = 1),
    error = function(e) e
  )
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "low acceptance rate")
  expect_no_match(conditionMessage(err), "cannot be started")
})

test_that("a non-finite or non-positive n_pool is refused, not left to size a batch wrongly", {
  # GitHub #190, same gap as max_steps above: n_pool sizes
  # `batch <- max(n_pool, n_chains)` in mcmc_init_resample()/mcmc_init_proposal()
  # with no validation.
  prior <- prior_uniform(low = -10, high = 10)
  lp <- function(theta) stats::dnorm(theta[, 1], mean = 5, sd = 0.5, log = TRUE)

  expect_error(
    mcmc_init(prior, lp, n_chains = 4, n_pool = NA),
    "`n_pool`"
  )
  expect_error(
    mcmc_init(prior, lp, n_chains = 4, n_pool = 0),
    "`n_pool`"
  )
  expect_error(
    mcmc_init(prior, lp, n_chains = 4, n_pool = -1),
    "`n_pool`"
  )
  expect_error(
    mcmc_init(prior, lp, n_chains = 4, strategy = "proposal", n_pool = NA),
    "`n_pool`"
  )
})

test_that("split-Rhat flags chains that disagree", {
  set.seed(6)
  converged <- array(stats::rnorm(400 * 4), c(400, 4, 1))
  # Same spread, different centers: exactly what Rhat is for.
  offset <- converged
  for (c in 1:4) offset[, c, 1] <- offset[, c, 1] + c

  expect_lt(mcmc_diagnostics(converged)$rhat, 1.01)
  expect_gt(mcmc_diagnostics(offset)$rhat, 1.5)
})

test_that("bulk ESS falls as autocorrelation rises", {
  set.seed(7)
  ar1 <- function(rho, n = 500, k = 4) {
    out <- array(0, c(n, k, 1))
    for (c in seq_len(k)) {
      y <- numeric(n)
      for (i in 2:n) y[i] <- rho * y[i - 1] + stats::rnorm(1, sd = sqrt(1 - rho^2))
      out[, c, 1] <- y
    }
    out
  }
  ess <- vapply(c(0.3, 0.5, 0.8, 0.95),
                function(r) mcmc_diagnostics(ar1(r))$ess_bulk, numeric(1))

  expect_true(all(diff(ess) < 0))
  # For an AR(1) chain tau = (1 + rho) / (1 - rho), so ESS ~ n * k / tau.
  expect_equal(ess, 2000 / ((1 + c(0.3, 0.5, 0.8, 0.95)) /
                              (1 - c(0.3, 0.5, 0.8, 0.95))),
               tolerance = 0.15)
})

test_that("rhat and bulk ESS agree with the posterior package", {
  skip_if_no_posterior()
  set.seed(8)
  n <- 400
  chains <- array(0, c(n, 4, 1))
  for (c in 1:4) {
    y <- numeric(n)
    for (i in 2:n) y[i] <- 0.7 * y[i - 1] + stats::rnorm(1)
    chains[, c, 1] <- y
  }
  ours <- mcmc_diagnostics(chains)

  # split_rhat() is max(bulk-Rhat, tail-Rhat) computed the same way
  # posterior::rhat() computes it by default, so this is the same
  # computation, not merely a close one: on this fixture the two agree to
  # floating-point precision.
  expect_equal(ours$rhat, posterior::rhat(chains[, , 1]), tolerance = 1e-6)
  expect_equal(ours$ess_bulk, posterior::ess_bulk(chains[, , 1]),
               tolerance = 0.01)
})

test_that("split_rhat() catches disagreement a classical Rhat misses", {
  # Rank-normalization is what makes Rhat robust to heavy tails: an extreme
  # outlier moves a raw value by an unbounded amount but a rank by at most
  # one, so a few huge draws cannot swamp a real between-chain difference the
  # way they can when Rhat is computed on raw values. Four chains centered at
  # 1, 2, 3, 4 with heavy Cauchy noise reproduce exactly that: the raw-value
  # (classical, pre-fix) statistic reads as converged because a handful of
  # huge draws inflate its within-chain variance enough to hide the offset,
  # while the rank-normalized one still sees the chains disagree.
  set.seed(42)
  n <- 400
  k <- 4
  m <- matrix(0, n, k)
  for (c in seq_len(k)) m[, c] <- c + stats::rcauchy(n, scale = 0.3)

  classical_rhat <- function(m) {
    chain_means <- colMeans(m)
    chain_vars <- apply(m, 2, stats::var)
    W <- mean(chain_vars)
    B <- nrow(m) * stats::var(chain_means)
    var_hat <- ((nrow(m) - 1) / nrow(m)) * W + B / nrow(m)
    sqrt(var_hat / W)
  }

  expect_lt(classical_rhat(m), 1.01)
  expect_gt(mcmc_diagnostics(array(m, c(n, k, 1)))$rhat, 1.2)
})

# GitHub #154 follow-up: the first fix landed with bulk-Rhat only, rank-
# normalizing the raw draws but never folding them, so split_rhat() still
# missed the failure mode tail-Rhat exists for -- chains that agree in
# location but disagree in spread. Bulk-Rhat alone reads such chains as
# converged because rank-normalizing the raw values still tracks location,
# not scale; folding around the median before rank-normalizing is what makes
# the statistic sensitive to a spread mismatch.
test_that("split_rhat() catches chains that agree in location but disagree in spread", {
  set.seed(154)
  n <- 400
  k <- 4
  m <- matrix(0, n, k)
  m[, 1] <- stats::rnorm(n, 0, 1)
  m[, 2] <- stats::rnorm(n, 0, 1)
  m[, 3] <- stats::rnorm(n, 0, 4)
  m[, 4] <- stats::rnorm(n, 0, 4)

  # What bulk-Rhat alone computes: rank-normalize the raw draws, no folding.
  bulk_only <- gelman_rubin_rhat(rank_normalize(m))
  expect_lt(bulk_only, 1.01)

  # split_rhat() adds the folded tail component and catches it.
  expect_gt(split_rhat(m), 1.2)
})

test_that("format_mcmc_diagnostics() reports a run it could not score", {
  # split_rhat() and bulk_ess() return NA for a run with too few iterations or
  # one chain, and max(na.rm = TRUE) over all-NA is -Inf with a warning. That
  # -Inf used to be printed as a diagnostic value.
  degenerate <- mcmc_diagnostics(array(stats::rnorm(6), dim = c(3L, 1L, 2L)))
  expect_true(all(is.na(degenerate$rhat)))
  expect_match(format_mcmc_diagnostics(degenerate), "diagnostics unavailable")
  expect_false(grepl("Inf", format_mcmc_diagnostics(degenerate), fixed = TRUE))

  # Some parameters scored and some not: report the ones that were, and count
  # the ones that were not.
  partial <- data.frame(rhat = c(1.01, NA_real_), ess_bulk = c(420, NA_real_))
  expect_match(format_mcmc_diagnostics(partial), "max Rhat 1.010")
  expect_match(format_mcmc_diagnostics(partial), "min bulk ESS 420")
  expect_match(format_mcmc_diagnostics(partial), "1 parameter not scored")

  # One statistic available and the other not.
  no_ess <- data.frame(rhat = c(1.2, 1.0), ess_bulk = c(NA_real_, NA_real_))
  expect_match(format_mcmc_diagnostics(no_ess), "max Rhat 1.200")
  expect_match(format_mcmc_diagnostics(no_ess), "bulk ESS unavailable")

  fine <- data.frame(rhat = c(1.0, 1.02), ess_bulk = c(900, 750))
  expect_identical(format_mcmc_diagnostics(fine),
                   "max Rhat 1.020, min bulk ESS 750")
})

test_that("split_rhat() returns NA for too few draws, too few chains, or a degenerate run", {
  # mcmc_diagnostics() never reaches these itself -- its own n_iter < 4 guard
  # and the always-even split-chain count keep split_rhat()'s shapes above
  # this floor -- but split_rhat() is `@keywords internal` and called
  # directly elsewhere, so it validates on its own.
  expect_true(is.na(split_rhat(matrix(1:4, nrow = 1))))  # n = 1
  expect_true(is.na(split_rhat(matrix(1:4, ncol = 1))))  # k = 1
  # a chain with a non-finite variance (here, an NA draw)
  expect_true(is.na(split_rhat(matrix(c(1, 2, NA, 4, 5, 6), nrow = 3, ncol = 2))))
  # every chain constant: zero within-chain variance
  expect_true(is.na(split_rhat(matrix(1, nrow = 5, ncol = 4))))
})

test_that("bulk_ess() returns NA for too few draws or a degenerate run", {
  expect_true(is.na(bulk_ess(matrix(1:8, nrow = 2, ncol = 4))))  # n = 2 < 4
  # every chain constant: rank-normalized variance is 0
  expect_true(is.na(bulk_ess(matrix(1, nrow = 5, ncol = 4))))
})
