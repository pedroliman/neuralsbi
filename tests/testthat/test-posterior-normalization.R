# Leakage correction: with a bounded prior the estimator places mass outside
# the support; log_prob(normalize = TRUE) must renormalize so the density
# integrates to one over the support and is -Inf outside it.

test_that("normalized log_prob integrates to one over a bounded support", {
  set.seed(11)
  # tight uniform prior so a nontrivial fraction of estimator mass leaks out
  prior <- prior_uniform(-1, 1)
  simulator <- function(theta) {
    theta + matrix(rnorm(length(theta), sd = 0.8), nrow = nrow(theta))
  }
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
  simulator <- function(theta) {
    theta + matrix(rnorm(length(theta), sd = 0.5), nrow = nrow(theta))
  }
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
  simulator <- function(theta) {
    theta + matrix(rnorm(length(theta), sd = 0.5), nrow = nrow(theta))
  }
  fit <- npe(prior, simulator, n_simulations = 1000,
             density_estimator = "linear_gaussian")
  post <- posterior(fit, x_obs = 0.5)
  theta <- matrix(seq(-2, 2, length.out = 11), ncol = 1)
  expect_equal(log_prob(post, theta, normalize = TRUE),
               log_prob(post, theta, normalize = FALSE))
})
