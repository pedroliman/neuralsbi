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
  # linear_gaussian, so this is the "Provide either" message from
  # prepare_simulations() and not check_torch_for_estimator()'s (#250) --
  # the two are independent checks and this test is about the former.
  expect_error(npe(prior, density_estimator = "linear_gaussian"), "simulator")
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

test_that("npe() errors on a flattened multi-column theta instead of guessing its layout (#291)", {
  # A column-major as.vector() flatten of an n x d theta matrix is a valid
  # multiple of d, so it slipped past both check_numeric() and the old
  # as_theta_matrix() call: matrix(theta_flat, ncol = d, byrow = TRUE)
  # silently reinterpreted it as a row-major layout, training on scrambled
  # parameter values with no warning at all.
  set.seed(1)
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  n <- 50
  theta_true <- matrix(c(rnorm(n, 0), rnorm(n, 5)), ncol = 2)
  theta_flat <- as.vector(theta_true)
  x <- theta_true + matrix(rnorm(n * 2, sd = 0.1), ncol = 2)

  expect_error(
    npe(prior, theta = theta_flat, x = x, density_estimator = "linear_gaussian"),
    "`theta` must be a matrix or data frame with 2 columns.*length-100 vector"
  )
  # nle() shares the same prepare_simulations() path.
  expect_error(
    nle(prior, theta = theta_flat, x = x, density_estimator = "linear_gaussian"),
    "`theta` must be a matrix or data frame with 2 columns"
  )
})

test_that("npe() still accepts a real theta matrix for a multi-parameter prior", {
  # The fix for #291 must not touch legitimate matrix input: an actual n x d
  # matrix carries no row-major/column-major ambiguity and should keep working.
  set.seed(2)
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  n <- 50
  theta <- matrix(c(rnorm(n, 0), rnorm(n, 5)), ncol = 2)
  x <- theta + matrix(rnorm(n * 2, sd = 0.1), ncol = 2)

  fit <- npe(prior, theta = theta, x = x, density_estimator = "linear_gaussian")
  expect_s3_class(fit, "nsbi_npe")
  expect_equal(unname(fit$std_theta$center), colMeans(theta), tolerance = 1e-8)
})

test_that("npe() still accepts a bare theta vector for a single-parameter prior", {
  # d == 1 has no row-major/column-major ambiguity (there is only one
  # column), so a bare vector must keep working exactly as before.
  set.seed(3)
  prior <- prior_normal(mean = 0, sd = 1)
  theta <- rnorm(200)
  x <- theta + rnorm(200, sd = 0.1)

  fit <- npe(prior, theta = theta, x = x, density_estimator = "linear_gaussian")
  expect_s3_class(fit, "nsbi_npe")
  expect_equal(fit$n_simulations, 200L)
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
