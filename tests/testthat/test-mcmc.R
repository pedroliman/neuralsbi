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

test_that("mcmc_init reports an impossible start rather than looping", {
  prior <- prior_uniform(low = 0, high = 1)
  expect_error(
    mcmc_init(prior, function(theta) rep(-Inf, nrow(theta)), n_chains = 2),
    "finite posterior density"
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

test_that("bulk ESS agrees with the posterior package", {
  skip_if_not_installed("posterior")
  set.seed(8)
  n <- 400
  chains <- array(0, c(n, 4, 1))
  for (c in 1:4) {
    y <- numeric(n)
    for (i in 2:n) y[i] <- 0.7 * y[i - 1] + stats::rnorm(1)
    chains[, c, 1] <- y
  }
  ours <- mcmc_diagnostics(chains)

  expect_equal(ours$rhat, posterior::rhat(chains[, , 1]), tolerance = 0.01)
  expect_equal(ours$ess_bulk, posterior::ess_bulk(chains[, , 1]),
               tolerance = 0.05)
})
