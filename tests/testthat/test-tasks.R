# Benchmark task definitions: shapes, priors, and analytic references.

test_that("gaussian_linear task matches its analytic reference", {
  set.seed(9)
  task <- task_gaussian_linear(dim = 3L)
  expect_s3_class(task, "nsbi_task")
  sims <- simulate_for_sbi(task$simulator, task$prior, 2000, seed = 9)
  expect_equal(dim(sims$x), c(2000L, 3L))

  # NPE with the exact estimator should match the analytic posterior.
  # Use an observation well away from zero so the posterior means are not
  # near-zero (where relative tolerance would be misleading).
  fit <- npe(task$prior, theta = sims$theta, x = sims$x,
             density_estimator = "linear_gaussian")
  x_obs <- c(0.3, -0.25, 0.2)
  draws <- sample(posterior(fit, x_obs = x_obs), 20000)
  ref <- task$reference_posterior(x_obs, 20000)
  # absolute differences: posterior sd here is ~0.22, so 0.03 is a tight bar
  expect_lt(max(abs(colMeans(draws) - colMeans(ref))), 0.03)
  expect_lt(max(abs(apply(draws, 2, sd) - apply(ref, 2, sd))), 0.03)
})

test_that("two_moons simulator produces the crescent geometry", {
  set.seed(10)
  task <- task_two_moons()
  theta <- matrix(0, nrow = 500, ncol = 2)  # 500 draws at fixed parameters
  x <- run_simulator(task$simulator, theta)
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
  x_fast <- task$simulator(c(0.9, 0.1))
  x_slow <- task$simulator(c(0.15, 0.14))
  expect_gt(max(x_fast), max(x_slow))
})

test_that("sir integrator keeps S and I within [0, N] for the prior's tail (#357)", {
  # Regression test for the Euler-overshoot bug: a single dt = 1 day step let
  # newinf exceed S for beta draws in the prior's tail, driving S negative
  # and I above N. That silently produced a physically-impossible trajectory
  # whose *output* still passed the pmin/pmax([0, 1]) clamp downstream, so
  # this checks the raw S/I path directly rather than the clamped output.
  # ~0.7-1% of draws hit the failure historically, so this draws enough theta
  # to reliably exercise it.
  set.seed(357)
  N <- 1e6; days <- 160
  prior <- prior_lognormal(meanlog = c(log(0.4), log(0.125)), sdlog = c(0.5, 0.2))
  theta <- sample_prior(prior, 2000)
  for (i in seq_len(nrow(theta))) {
    path <- sir_trajectory(theta[i, ], N = N, days = days)
    expect_true(all(path$S >= 0 & path$S <= N),
                info = sprintf("row %d: beta=%.4f gamma=%.4f", i, theta[i, 1], theta[i, 2]))
    expect_true(all(path$I >= 0 & path$I <= N),
                info = sprintf("row %d: beta=%.4f gamma=%.4f", i, theta[i, 1], theta[i, 2]))
  }
})

test_that("print.nsbi_task() reports dimensions and whether a reference posterior exists", {
  expect_output(print(task_gaussian_linear(dim = 3L)),
               "gaussian_linear: 3 parameters -> 3 data dims \\(analytic reference available\\)")
  expect_output(print(task_two_moons()),
               "two_moons: 2 parameters -> 2 data dims$")
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
  x <- run_simulator(task$simulator, theta_fix)
  expect_equal(mean(x[, 1]), 1, tolerance = 0.05)
  expect_equal(mean(x[, 4]), -1, tolerance = 0.05)
})
