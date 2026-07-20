# TSNPE with the exact linear_gaussian estimator: the truncation machinery is
# exercised end to end, and the final posterior must still match the analytic
# conjugate posterior at the targeted observation.

analytic_gauss_posterior <- function(x_obs, sigma, d) {
  prec <- diag(d) + diag(d) / sigma^2
  Sigma <- solve(prec)
  mu <- as.numeric(Sigma %*% (x_obs / sigma^2))
  list(mu = mu, Sigma = Sigma)
}

test_that("npe_sequential recovers the analytic posterior at x_obs", {
  set.seed(21)
  d <- 2; sigma <- 0.5
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) {
    theta + matrix(rnorm(length(theta), sd = sigma), nrow = nrow(theta))
  }
  x_obs <- c(1.0, -0.5)
  fit <- npe_sequential(prior, simulator, x_obs = x_obs,
                        n_rounds = 3, n_simulations = 800,
                        density_estimator = "linear_gaussian", seed = 1)
  expect_s3_class(fit, "nsbi_snpe")
  expect_s3_class(fit, "nsbi_npe")
  expect_length(fit$rounds, 3L)
  expect_equal(fit$n_simulations, 2400L)

  post <- posterior(fit, x_obs = x_obs)
  draws <- sample(post, 5000)
  truth <- analytic_gauss_posterior(x_obs, sigma, d)
  expect_equal(colMeans(draws), truth$mu, tolerance = 0.06)
  expect_equal(apply(draws, 2, sd), sqrt(diag(truth$Sigma)), tolerance = 0.06)

  z <- matrix(rnorm(5000 * d), ncol = d)
  analytic_draws <- sweep(z %*% chol(truth$Sigma), 2, truth$mu, `+`)
  expect_lt(c2st(draws, analytic_draws, seed = 2)$accuracy, 0.6)
})

test_that("later rounds actually truncate the proposal", {
  set.seed(22)
  prior <- prior_normal(mean = 0, sd = 2)
  simulator <- function(theta) {
    theta + matrix(rnorm(length(theta), sd = 0.2), nrow = nrow(theta))
  }
  fit <- npe_sequential(prior, simulator, x_obs = 1,
                        n_rounds = 2, n_simulations = 500,
                        density_estimator = "linear_gaussian", seed = 3)
  # round 1 draws from the prior itself
  expect_equal(fit$rounds[[1]]$acceptance, 1)
  expect_identical(fit$rounds[[1]]$threshold, -Inf)
  # posterior sd ~0.2 vs prior sd 2: most prior candidates must be rejected
  expect_lt(fit$rounds[[2]]$acceptance, 0.8)
  expect_gt(fit$rounds[[2]]$acceptance, 0)
  expect_true(is.finite(fit$rounds[[2]]$threshold))
  expect_output(print(fit), "nsbi_snpe")
})

test_that("per-round budgets can differ and are recorded", {
  set.seed(23)
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) {
    theta + matrix(rnorm(length(theta), sd = 0.5), nrow = nrow(theta))
  }
  fit <- npe_sequential(prior, simulator, x_obs = 0.3,
                        n_rounds = 2, n_simulations = c(600, 300),
                        density_estimator = "linear_gaussian", seed = 4)
  expect_equal(fit$rounds[[1]]$n_new, 600L)
  expect_equal(fit$rounds[[2]]$n_new, 300L)
  expect_equal(fit$n_simulations, 900L)
})

test_that("npe_sequential requires a simulator function", {
  prior <- prior_normal(0, 1)
  expect_error(npe_sequential(prior, simulator = NULL, x_obs = 0),
               "simulator")
})

test_that("truncation works with a bounded prior", {
  set.seed(24)
  prior <- prior_uniform(c(-2, -2), c(2, 2))
  simulator <- function(theta) {
    theta + matrix(rnorm(length(theta), sd = 0.3), nrow = nrow(theta))
  }
  x_obs <- c(0.5, 0.5)
  fit <- npe_sequential(prior, simulator, x_obs = x_obs,
                        n_rounds = 2, n_simulations = 600,
                        density_estimator = "linear_gaussian", seed = 5)
  post <- posterior(fit, x_obs = x_obs)
  draws <- sample(post, 2000)
  # posterior concentrated near the truth, well inside the box
  expect_equal(colMeans(draws), x_obs, tolerance = 0.1)
  expect_true(all(within_support(prior, draws)))
})
