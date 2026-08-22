# The linear_gaussian estimator is exact for a linear-Gaussian model, so it
# gives a torch-free regression oracle for the whole NPE pipeline.

test_that("linear_gaussian NPE recovers the analytic Gaussian posterior", {
  set.seed(42)
  d <- 2; sigma <- 0.5
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + rnorm(length(theta), sd = sigma)
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
  expect_lt(c2st(analytic_draws, draws, classifier = "logistic",
                 seed = 1)$accuracy, 0.6)
  # The sbibm classifier also sees a difference in spread, which the linear
  # one cannot. It needs torch, so it runs only where libtorch is installed.
  if (has_torch()) {
    expect_lt(c2st(analytic_draws, draws, seed = 1)$accuracy, 0.6)
  }
})

test_that("npe errors clearly when neither simulator nor (theta,x) are given", {
  prior <- prior_normal(0, 1)
  expect_error(npe(prior), "simulator")
})

test_that("a bad density_estimator errors before the simulator runs", {
  prior <- prior_normal(mean = 0, sd = 1)
  calls <- 0L
  counting_simulator <- function(theta) {
    calls <<- calls + 1L
    theta
  }
  expect_error(
    npe(prior, counting_simulator, n_simulations = 100,
        density_estimator = "mfa"),
    "should be one of"
  )
  expect_identical(calls, 0L)
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

test_that("a non-finite theta is dropped on the pre-computed path", {
  set.seed(11)
  prior <- prior_normal(mean = 0, sd = 1)
  theta <- sample_prior(prior, 200)
  x <- theta + matrix(rnorm(200, sd = 0.3), ncol = 1)
  theta[3, 1] <- NA

  # without this the NA reaches chol() in the estimator and comes back as
  # "the leading minor of order 1 is not positive"
  expect_warning(
    fit <- npe(prior, theta = theta, x = x,
               density_estimator = "linear_gaussian"),
    "Dropped 1 of 200 simulations with non-finite parameters"
  )
  expect_equal(fit$n_simulations, 199L)
  expect_equal(fit$n_dropped, 1L)

  # a bad theta and a bad x on different rows both go
  x[7, 1] <- Inf
  expect_warning(
    fit2 <- npe(prior, theta = theta, x = x,
                density_estimator = "linear_gaussian"),
    "Dropped 2 of 200 simulations with non-finite parameters or output"
  )
  expect_equal(fit2$n_simulations, 198L)

  # nothing left is an error, and it points at theta rather than the simulator
  x[7, 1] <- 0
  theta[] <- NA_real_
  expect_error(
    npe(prior, theta = theta, x = x, density_estimator = "linear_gaussian"),
    "All 200 simulations returned non-finite parameters"
  )
})
