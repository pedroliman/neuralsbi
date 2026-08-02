# An NPE posterior conditions on one observation, so resolve_x() keeps row 1.
# The rows of an NLE observation are independent draws instead, so the same
# multi-row `x_obs` means something different to each. The warning is what
# tells a user which of the two they got.

npe_fit_2d <- function(seed = 1) {
  set.seed(seed)
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + stats::rnorm(2, sd = 0.5)
  npe(prior, simulator, n_simulations = 400,
      density_estimator = "linear_gaussian")
}

test_that("a multi-row observation warns and names the row count", {
  fit <- npe_fit_2d()
  x_obs <- matrix(stats::rnorm(40), ncol = 2)
  post <- posterior(fit, x_obs = x_obs)

  expect_warning(sample(post, 20), "20 rows")
  expect_warning(sample(post, 20), "row 1 is used")
  expect_warning(sample(post, 20), "nle\\(\\)")
  expect_warning(log_prob(post, c(0, 0)), "20 rows")
  # also when the observation comes in through the call rather than the object
  post_none <- posterior(fit)
  expect_warning(sample(post_none, 20, obs = x_obs), "20 rows")
})

test_that("the warning does not change which row is used", {
  fit <- npe_fit_2d()
  x_obs <- matrix(stats::rnorm(40), ncol = 2)
  post_multi <- posterior(fit, x_obs = x_obs)
  post_one <- posterior(fit, x_obs = x_obs[1, , drop = FALSE])

  lp_multi <- suppressWarnings(log_prob(post_multi, c(0.2, -0.3)))
  expect_equal(lp_multi, log_prob(post_one, c(0.2, -0.3)))
})

test_that("a single row or a plain observation vector is silent", {
  fit <- npe_fit_2d()
  post_vec <- posterior(fit, x_obs = c(0.4, -0.4))
  post_mat <- posterior(fit, x_obs = matrix(c(0.4, -0.4), nrow = 1))

  expect_silent(sample(post_vec, 20))
  expect_silent(sample(post_mat, 20))
  expect_silent(log_prob(post_vec, c(0, 0)))
  expect_silent(log_prob(post_mat, c(0, 0)))
})

test_that("an NLE posterior does not warn about repeated observations", {
  set.seed(2)
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) stats::rnorm(1, theta, sd = 0.5)
  fit <- nle(prior, simulator, n_simulations = 500,
             density_estimator = "linear_gaussian", seed = 2)
  x_obs <- matrix(stats::rnorm(20, mean = 0.3, sd = 0.5), ncol = 1)
  post <- posterior(fit, x_obs, n_chains = 4, warmup = 20, thin = 1, seed = 3)

  expect_no_warning(sample(post, 50))
  expect_no_warning(log_prob(post, 0.3))
})
