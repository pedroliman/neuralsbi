# Leakage correction: with a bounded prior the estimator places mass outside
# the support; log_prob(normalize = TRUE) must renormalize so the density
# integrates to one over the support and is -Inf outside it.

test_that("normalized log_prob integrates to one over a bounded support", {
  set.seed(11)
  # tight uniform prior so a nontrivial fraction of estimator mass leaks out
  prior <- prior_uniform(-1, 1)
  simulator <- function(theta) theta + rnorm(length(theta), sd = 0.8)
  fit <- npe(prior, simulator, n_simulations = 3000,
             density_estimator = "linear_gaussian")
  post <- posterior(fit, x_obs = 0.9)

  grid <- matrix(seq(-1, 1, length.out = 2001), ncol = 1)
  h <- diff(grid[1:2, 1])

  lp_raw <- log_prob(post, grid, normalize = FALSE)
  mass_raw <- sum(exp(lp_raw)) * h
  # the unnormalized density leaks mass outside [-1, 1]
  expect_lt(mass_raw, 0.99)

  lp <- log_prob(post, grid, normalize = TRUE)
  mass <- sum(exp(lp)) * h
  expect_equal(mass, 1, tolerance = 0.05)
})

test_that("normalized log_prob is -Inf outside the support and shifted inside", {
  set.seed(12)
  prior <- prior_uniform(c(-1, -1), c(1, 1))
  simulator <- function(theta) theta + rnorm(length(theta), sd = 0.5)
  fit <- npe(prior, simulator, n_simulations = 2000,
             density_estimator = "linear_gaussian")
  post <- posterior(fit, x_obs = c(0.8, 0.8))

  outside <- rbind(c(2, 0), c(0, -1.5), c(5, 5))
  expect_true(all(log_prob(post, outside) == -Inf))
  # unnormalized evaluation still returns finite densities out there
  expect_true(all(is.finite(log_prob(post, outside, normalize = FALSE))))

  inside <- rbind(c(0.5, 0.5), c(0, 0), c(-0.9, 0.9))
  lp_raw <- log_prob(post, inside, normalize = FALSE)
  lp <- log_prob(post, inside)
  shift <- lp - lp_raw
  # renormalization adds the same -log(acceptance) >= 0 everywhere inside
  expect_true(all(shift >= 0))
  expect_equal(max(shift) - min(shift), 0, tolerance = 1e-12)
})

test_that("unbounded priors are unaffected by normalize", {
  set.seed(13)
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + rnorm(length(theta), sd = 0.5)
  fit <- npe(prior, simulator, n_simulations = 1000,
             density_estimator = "linear_gaussian")
  post <- posterior(fit, x_obs = 0.5)
  theta <- matrix(seq(-2, 2, length.out = 11), ncol = 1)
  expect_equal(log_prob(post, theta, normalize = TRUE),
               log_prob(post, theta, normalize = FALSE))
})

test_that("sample() requests the shortfall each rejection-sampling round, not a full batch, and still returns exactly n rows", {
  set.seed(11)
  prior <- prior_uniform(low = c(0, 0), high = c(1, 1))
  simulator <- function(theta) theta + rnorm(length(theta), sd = 0.2)
  fit <- npe(prior, simulator, n_simulations = 500,
             density_estimator = "linear_gaussian")
  # shift the fitted mean partway out of the box: leaks enough mass that more
  # than one rejection-sampling round is needed, without exhausting
  # max_sampling_batches
  fit$de$B[1, ] <- fit$de$B[1, ] + 0.5
  post <- posterior(fit, x_obs = c(0.9, 0.9))

  requested <- integer(0)
  real_de_sample <- de_sample
  local_mocked_bindings(
    de_sample = function(de, x, n) {
      requested[[length(requested) + 1L]] <<- n
      real_de_sample(de, x, n)
    }
  )

  draws <- sample(post, n = 500, max_sampling_batches = 50)
  expect_equal(nrow(draws), 500L)
  expect_gt(length(requested), 1L)   # more than one round was needed here
  expect_equal(requested[1], 500L)   # first round always asks for the full n
  # every later round asks only for the still-missing draws, never the full
  # 500 again, and the shortfall never grows round to round
  expect_true(all(requested[-1] < 500L))
  expect_true(all(diff(requested) <= 0))
})

test_that("the posterior counts are checked before they reach de_sample()", {
  set.seed(12)
  prior <- prior_uniform(-1, 1)
  simulator <- function(theta) theta + stats::rnorm(1, sd = 0.5)
  fit <- npe(prior, simulator, n_simulations = 300,
             density_estimator = "linear_gaussian")
  post <- posterior(fit, x_obs = 0.2)

  expect_error(sample(post, n = 100, max_sampling_batches = 0),
               "`max_sampling_batches` must be a single whole number of at least 1 since")
  expect_error(sample(post, n = 100, max_sampling_batches = 2.5),
               "not 2.5")
  expect_error(log_prob(post, 0.2, n_normalization = 0),
               "`n_normalization` must be")
  expect_error(log_prob(post, 0.2, n_normalization = NA),
               "`n_normalization` must be")
  expect_error(map_estimate(post, n_init = 0),
               "`n_init` must be a single whole number of at least 1 since")

  # An unbounded prior never reads n_normalization, but the argument is
  # checked either way: a wrong value is a mistake in the call whichever
  # branch this run takes.
  unbounded <- npe(prior_normal(mean = 0, sd = 1), simulator,
                   n_simulations = 300, density_estimator = "linear_gaussian")
  expect_error(log_prob(posterior(unbounded, x_obs = 0.2), 0.2,
                        n_normalization = -1),
               "`n_normalization` must be")
})

# GitHub #152: map_estimate() optimizes with unconstrained Nelder-Mead and
# always calls log_prob(normalize = FALSE), so a bounded prior's -Inf mask
# never reached the objective and the search could walk outside the box.
test_that("map_estimate() stays inside a bounded prior's support", {
  set.seed(14)
  prior <- prior_uniform(0, 1)
  simulator <- function(theta) theta + stats::rnorm(1, sd = 0.05)
  fit <- npe(prior, simulator, n_simulations = 500,
             density_estimator = "linear_gaussian")
  # nudge the fitted conditional mean outside [0, 1] so the unconstrained
  # optimum sits outside the prior's support
  fit$de$B[1, ] <- fit$de$B[1, ] + 0.5
  post <- posterior(fit, x_obs = 0.9)

  map <- suppressWarnings(map_estimate(post))
  expect_true(within_support(prior, matrix(map, nrow = 1)))
})

test_that("map_estimate() respects a one-sided (lower-only) bound", {
  set.seed(15)
  prior <- prior_custom(
    sample_fn = function(n) matrix(stats::rexp(n, 1), ncol = 1),
    log_prob_fn = function(theta) stats::dexp(theta[, 1], 1, log = TRUE),
    dim = 1, lower = 0
  )
  simulator <- function(theta) theta + stats::rnorm(1, sd = 0.3)
  fit <- npe(prior, simulator, n_simulations = 2000,
             density_estimator = "linear_gaussian")
  # nudge the fitted conditional mean below the lower bound
  fit$de$B[1, ] <- fit$de$B[1, ] - 1
  post <- posterior(fit, x_obs = 0.1)

  map <- suppressWarnings(map_estimate(post))
  expect_true(within_support(prior, matrix(map, nrow = 1)))
})

test_that("map_estimate() on an unbounded prior lands between the prior mean and the observation", {
  set.seed(13)
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + stats::rnorm(1, sd = 0.5)
  fit <- npe(prior, simulator, n_simulations = 1000,
             density_estimator = "linear_gaussian")
  post <- posterior(fit, x_obs = 0.5)

  map <- suppressWarnings(map_estimate(post))
  expect_true(is.finite(map))
  expect_true(map > 0 && map < 0.5)
})
