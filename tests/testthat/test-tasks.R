# Benchmark task definitions: shapes, priors, and analytic references.

test_that("gaussian_linear task matches its analytic reference", {
  set.seed(9)
  task <- task_gaussian_linear(dim = 3L)
  expect_s3_class(task, "nsbi_task")
  sims <- simulate_for_sbi(task$simulator, task$prior, 2000, seed = 9)
  expect_equal(dim(sims$x), c(2000L, 3L))

  # NPE with the exact estimator should match the analytic posterior
  fit <- npe(task$prior, theta = sims$theta, x = sims$x,
             density_estimator = "linear_gaussian")
  x_obs <- c(0.1, -0.2, 0.05)
  draws <- sample(posterior(fit, x_obs = x_obs), 5000)
  ref <- task$reference_posterior(x_obs, 5000)
  expect_equal(colMeans(draws), colMeans(ref), tolerance = 0.05)
  expect_equal(apply(draws, 2, sd), apply(ref, 2, sd), tolerance = 0.05)
})

test_that("two_moons simulator produces the crescent geometry", {
  set.seed(10)
  task <- task_two_moons()
  theta <- matrix(0, nrow = 500, ncol = 2)  # fixed parameters
  x <- task$simulator(theta)
  expect_equal(dim(x), c(500L, 2L))
  # crescent radius approx 0.1 around (0.25, 0)
  r <- sqrt((x[, 1] - 0.25)^2 + x[, 2]^2)
  expect_equal(mean(r), 0.1, tolerance = 0.02)
  expect_true(all(x[, 1] >= 0.25 - 0.2))  # right half-moon only
})

test_that("sir task simulates plausible epidemics under its prior", {
  set.seed(12)
  task <- task_sir()
  sims <- simulate_for_sbi(task$simulator, task$prior, 50)
  expect_equal(dim(sims$x), c(50L, 10L))
  expect_true(all(sims$x >= 0 & sims$x <= 1))       # infected fractions
  expect_true(all(within_support(task$prior, sims$theta)))
  expect_true(all(sims$theta > 0))                   # rates are positive
  lp <- task$prior$log_prob(sims$theta)
  expect_true(all(is.finite(lp)))
  # a fast-spreading epidemic peaks higher than a slow one
  x_fast <- task$simulator(matrix(c(0.9, 0.1), nrow = 1))
  x_slow <- task$simulator(matrix(c(0.15, 0.14), nrow = 1))
  expect_gt(max(x_fast), max(x_slow))
})

test_that("slcp simulator has the right shape and support behavior", {
  set.seed(11)
  task <- task_slcp()
  sims <- simulate_for_sbi(task$simulator, task$prior, 100)
  expect_equal(dim(sims$theta), c(100L, 5L))
  expect_equal(dim(sims$x), c(100L, 8L))
  expect_true(all(within_support(task$prior, sims$theta)))
  # the four 2-D points are i.i.d.: odd columns share mean theta_1
  theta_fix <- matrix(rep(c(1, -1, 0.8, 0.6, 0.3), each = 2000), ncol = 5)
  x <- task$simulator(theta_fix)
  expect_equal(mean(x[, 1]), 1, tolerance = 0.05)
  expect_equal(mean(x[, 4]), -1, tolerance = 0.05)
})
