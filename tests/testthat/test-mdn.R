# These tests exercise the neural (torch) MDN. They are skipped automatically
# when torch / libtorch is not installed (see helper-torch.R), so the suite
# still runs everywhere.

test_that("MDN log_prob and sampling shapes are correct", {
  skip_if_no_torch()
  set.seed(1)
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + rnorm(length(theta), sd = 0.5)
  fit <- npe(prior, simulator, n_simulations = 2000, density_estimator = "mdn",
             n_components = 2L, hidden = c(30L, 30L), max_epochs = 60L,
             seed = 1)
  post <- posterior(fit, x_obs = c(1, -0.5))
  draws <- sample(post, 500)
  expect_equal(dim(draws), c(500L, 2L))
  lp <- log_prob(post, draws[1:10, ], normalize = FALSE)
  expect_length(lp, 10L)
  expect_true(all(is.finite(lp)))
})

test_that("MDN sampling works when dim_theta == 1 and n_components > 1 (#211)", {
  # Regression test: de_sample.nsbi_de_mdn used to index the batch dim of
  # torch tensors with `[1, , ]`, which follows R's drop=TRUE semantics and
  # drops every size-1 dimension, not just the batch dim. With a 1-d target
  # and more than one mixture component, means/Larr collapsed to bare
  # vectors and de_sample() errored with "incorrect number of dimensions".
  skip_if_no_torch()
  set.seed(3)
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + rnorm(length(theta), sd = 0.5)
  fit <- npe(prior, simulator, n_simulations = 500, density_estimator = "mdn",
             n_components = 3L, hidden = c(16L, 16L), max_epochs = 10L,
             seed = 3)
  post <- posterior(fit, x_obs = 0.5)
  draws <- sample(post, 50)
  expect_equal(dim(draws), c(50L, 1L))
  expect_true(all(is.finite(draws)))
})

test_that("MDN posterior is close to the analytic linear-Gaussian posterior", {
  skip_if_no_torch()
  set.seed(2)
  d <- 2; sigma <- 0.5
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + rnorm(length(theta), sd = sigma)
  fit <- npe(prior, simulator, n_simulations = 8000, density_estimator = "mdn",
             n_components = 3L, hidden = c(50L, 50L), max_epochs = 300L, seed = 2)
  x_obs <- c(1.0, -0.5)
  post <- posterior(fit, x_obs = x_obs)
  draws <- sample(post, 10000)

  truth <- analytic_gauss_posterior(x_obs, sigma, d)
  mu <- truth$mu; Sigma <- truth$Sigma
  # neural fit: looser tolerance than the exact linear estimator
  expect_equal(colMeans(draws), mu, tolerance = 0.1)
  expect_equal(apply(draws, 2, sd), sqrt(diag(Sigma)), tolerance = 0.1)
})
