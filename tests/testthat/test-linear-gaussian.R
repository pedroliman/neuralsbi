# The linear_gaussian estimator is exact for a linear-Gaussian model, so it
# gives a torch-free regression oracle for the whole NPE pipeline.

analytic_gauss_posterior <- function(x_obs, sigma, d) {
  prec <- diag(d) + diag(d) / sigma^2
  Sigma <- solve(prec)
  mu <- as.numeric(Sigma %*% (x_obs / sigma^2))
  list(mu = mu, Sigma = Sigma)
}

test_that("linear_gaussian NPE recovers the analytic Gaussian posterior", {
  set.seed(42)
  d <- 2; sigma <- 0.5
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) {
    theta + matrix(rnorm(length(theta), sd = sigma), nrow = nrow(theta))
  }
  fit <- npe(prior, simulator, n_simulations = 4000,
             density_estimator = "linear_gaussian")
  x_obs <- c(1.0, -0.5)
  post <- posterior(fit, x_obs = x_obs)
  draws <- sample(post, 10000)

  truth <- analytic_gauss_posterior(x_obs, sigma, d)
  expect_equal(colMeans(draws), truth$mu, tolerance = 0.05)
  expect_equal(apply(draws, 2, sd), sqrt(diag(truth$Sigma)), tolerance = 0.05)

  # indistinguishable from analytic draws
  z <- matrix(rnorm(10000 * d), ncol = d)
  analytic_draws <- sweep(z %*% chol(truth$Sigma), 2, truth$mu, `+`)
  acc <- c2st(draws, analytic_draws, seed = 1)$accuracy
  expect_lt(acc, 0.6)
})

test_that("npe errors clearly when neither simulator nor (theta,x) are given", {
  prior <- prior_normal(0, 1)
  expect_error(npe(prior), "simulator")
})

test_that("pre-computed simulations can be passed directly", {
  prior <- prior_normal(mean = 0, sd = 1)
  theta <- sample_prior(prior, 500)
  x <- theta + matrix(rnorm(500, sd = 0.3), ncol = 1)
  fit <- npe(prior, theta = theta, x = x,
             density_estimator = "linear_gaussian")
  expect_s3_class(fit, "nsbi_npe")
  expect_equal(fit$n_simulations, 500L)
})
