test_that("sbc returns ranks of the right shape and reasonable calibration", {
  set.seed(7)
  d <- 2; sigma <- 0.5
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) {
    theta + matrix(rnorm(length(theta), sd = sigma), nrow = nrow(theta))
  }
  fit <- npe(prior, simulator, n_simulations = 4000,
             density_estimator = "linear_gaussian")
  res <- sbc(fit, simulator, n_sbc = 100, n_posterior_samples = 200, seed = 3)
  expect_equal(dim(res$ranks), c(100L, 2L))
  # a well-specified exact estimator should not fail uniformity badly
  expect_true(all(res$uniformity_pvalue > 0.001))
})

test_that("expected_coverage produces a monotone-ish curve near the diagonal", {
  set.seed(7)
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + matrix(rnorm(length(theta), sd = 0.5),
                                              nrow = nrow(theta))
  fit <- npe(prior, simulator, n_simulations = 4000,
             density_estimator = "linear_gaussian")
  res <- sbc(fit, simulator, n_sbc = 200, n_posterior_samples = 200, seed = 5)
  cov <- expected_coverage(res, levels = c(0.5, 0.9))
  expect_true(all(cov$param1 >= 0 & cov$param1 <= 1))
  # 90% interval should cover clearly more often than the 50% interval
  expect_gt(cov$param1[2], cov$param1[1])
})

test_that("c2st of a sample set against itself is ~0.5", {
  set.seed(1)
  a <- matrix(rnorm(2000), ncol = 2)
  b <- matrix(rnorm(2000), ncol = 2)
  res <- c2st(a, b, seed = 1)
  expect_lt(res$accuracy, 0.6)
})

test_that("posterior_predictive returns simulator-shaped output", {
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + matrix(rnorm(nrow(theta), sd = 0.3),
                                             ncol = 1)
  fit <- npe(prior, simulator, n_simulations = 1000,
             density_estimator = "linear_gaussian")
  post <- posterior(fit, x_obs = 0.5)
  pp <- posterior_predictive(post, simulator, n = 300)
  expect_equal(nrow(pp), 300L)
})

test_that("plot_posterior_predictive runs and locates the observation", {
  set.seed(2)
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + matrix(rnorm(length(theta), sd = 0.3),
                                              nrow = nrow(theta))
  fit <- npe(prior, simulator, n_simulations = 1000,
             density_estimator = "linear_gaussian")
  x_obs <- c(0.5, -0.2)
  post <- posterior(fit, x_obs = x_obs)
  pp <- posterior_predictive(post, simulator, n = 500)
  path <- tempfile(fileext = ".png")
  grDevices::png(path)
  q <- plot_posterior_predictive(pp, x_obs)
  grDevices::dev.off()
  expect_true(file.exists(path))
  expect_length(q, 2L)
  # the observation should sit inside the bulk of its own predictive
  expect_true(all(q > 0.01 & q < 0.99))
  unlink(path)
})
