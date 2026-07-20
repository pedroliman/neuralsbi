# These tests exercise the neural (torch) MDN. They are skipped automatically
# when torch / libtorch is not installed (see helper-torch.R), so the suite
# still runs everywhere.

test_that("MDN log_prob and sampling shapes are correct", {
  skip_if_no_torch()
  set.seed(1)
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + matrix(rnorm(length(theta), sd = 0.5),
                                             nrow = nrow(theta))
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

test_that("MDN posterior is close to the analytic linear-Gaussian posterior", {
  skip_if_no_torch()
  set.seed(2)
  d <- 2; sigma <- 0.5
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + matrix(rnorm(length(theta), sd = sigma),
                                             nrow = nrow(theta))
  fit <- npe(prior, simulator, n_simulations = 8000, density_estimator = "mdn",
             n_components = 3L, hidden = c(50L, 50L), max_epochs = 300L, seed = 2)
  x_obs <- c(1.0, -0.5)
  post <- posterior(fit, x_obs = x_obs)
  draws <- sample(post, 10000)

  prec <- diag(d) + diag(d) / sigma^2
  Sigma <- solve(prec)
  mu <- as.numeric(Sigma %*% (x_obs / sigma^2))
  # neural fit: looser tolerance than the exact linear estimator
  expect_equal(colMeans(draws), mu, tolerance = 0.1)
  expect_equal(apply(draws, 2, sd), sqrt(diag(Sigma)), tolerance = 0.1)
})
