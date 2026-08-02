# An observation is the one input a user types out rather than simulates, so it
# is the one most likely to carry an NA from a real data set. Left alone the NA
# runs through standardization into the estimator and comes back as all-NaN
# draws; the first complaint then comes from stats::quantile() inside summary(),
# which says nothing about the observation.

npe_fit_2d <- function(seed = 1) {
  set.seed(seed)
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + stats::rnorm(2, sd = 0.5)
  npe(prior, simulator, n_simulations = 400,
      density_estimator = "linear_gaussian")
}

nle_fit_1d <- function(seed = 2) {
  set.seed(seed)
  nle(prior_normal(mean = 0, sd = 1),
      function(theta) stats::rnorm(1, theta, sd = 0.5),
      n_simulations = 400, density_estimator = "linear_gaussian", seed = seed)
}

test_that("posterior() rejects a non-finite x_obs and names the entry", {
  fit <- npe_fit_2d()

  expect_error(posterior(fit, x_obs = c(0.5, NA)),
               "`x_obs` contains 1 non-finite value \\(NA\\), first at position 2")
  expect_error(posterior(fit, x_obs = c(NaN, 0.5)), "`x_obs` contains")
  expect_error(posterior(fit, x_obs = c(Inf, 0.5)),
               "`x_obs` contains 1 non-finite value \\(Inf\\)")
  expect_error(posterior(fit, x_obs = matrix(c(1, 2, 3, NA), ncol = 2)),
               "first at row 2, column 2")

  # a finite observation is still accepted, in either shape
  expect_s3_class(posterior(fit, x_obs = c(0.4, -0.4)), "nsbi_posterior")
  expect_s3_class(posterior(fit, x_obs = matrix(c(0.4, -0.4), nrow = 1)),
                  "nsbi_posterior")
})

test_that("an observation passed to sample() or log_prob() is checked too", {
  fit <- npe_fit_2d()
  post <- posterior(fit)

  expect_error(sample(post, 10, obs = c(0.5, NA)),
               "`x` contains 1 non-finite value \\(NA\\)")
  expect_error(log_prob(post, c(0, 0), x = c(NA, 0.5)),
               "`x` contains 1 non-finite value \\(NA\\)")
  expect_error(map_estimate(post, x = c(0.5, Inf)), "`x` contains")
})

test_that("an NLE posterior rejects a non-finite observation", {
  fit <- nle_fit_1d()

  expect_error(posterior(fit, x_obs = c(0.2, NA, 0.4)),
               "`x_obs` contains 1 non-finite value \\(NA\\), first at position 2")

  post <- posterior(fit, x_obs = matrix(stats::rnorm(10), ncol = 1),
                    n_chains = 4, warmup = 20, thin = 1, seed = 3)
  expect_error(sample(post, 20, obs = c(0.1, NA)),
               "`obs` contains 1 non-finite value \\(NA\\)")
  expect_error(log_prob(post, 0.3, x = c(NA, 0.1), normalize = FALSE),
               "`x` contains 1 non-finite value \\(NA\\)")
})
