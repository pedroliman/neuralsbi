test_that("sbc returns ranks of the right shape and reasonable calibration", {
  set.seed(7)
  d <- 2; sigma <- 0.5
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + rnorm(length(theta), sd = sigma)
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
  simulator <- function(theta) theta + rnorm(length(theta), sd = 0.5)
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

test_that("c2st is not fooled by unequal sample sizes", {
  # Accuracy against unbalanced classes is not a two-sample test. Four times as
  # many draws on one side and a classifier scores 0.8 by always answering with
  # the bigger class, having learned nothing -- which is what the vignette hit
  # comparing 8000 slice draws against 2000 from Stan.
  set.seed(2)
  a <- matrix(rnorm(16000), ncol = 2)
  b <- matrix(rnorm(4000), ncol = 2)

  expect_lt(c2st(a, b, seed = 1)$accuracy, 0.6)
  expect_lt(c2st(b, a, seed = 1)$accuracy, 0.6)

  # Balancing must not cost it the ability to see a real difference.
  shifted <- matrix(rnorm(4000, mean = 3), ncol = 2)
  expect_gt(c2st(a, shifted, seed = 1)$accuracy, 0.9)
})

test_that("posterior_predictive returns simulator-shaped output", {
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + rnorm(1, sd = 0.3)
  fit <- npe(prior, simulator, n_simulations = 1000,
             density_estimator = "linear_gaussian")
  post <- posterior(fit, x_obs = 0.5)
  pp <- posterior_predictive(post, simulator, n = 300)
  expect_equal(nrow(pp), 300L)
})

test_that("plot_posterior_predictive runs and locates the observation", {
  skip_if_no_ggplot2()
  set.seed(2)
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + rnorm(length(theta), sd = 0.3)
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
