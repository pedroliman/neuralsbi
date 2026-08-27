# GitHub #234: within_support() returns NA (not FALSE) for a row containing
# a NaN/NA value. log_prob.nsbi_posterior() already guards against this for
# user-supplied theta by running check_finite() first (#221, see
# test-posterior-nonfinite-theta.R), but three call sites test a vector that
# comes from the density estimator's own output -- de_sample() -- and had no
# such guard: sample.nsbi_posterior()'s rejection filter, log_prob()'s
# normalize = TRUE acceptance estimate, and map_estimate()'s objective. These
# tests mock de_sample() (or, for map_estimate(), stats::optim()) to inject a
# NaN draw the way an under-trained MAF/NSF/MDN might, and check that a
# non-finite draw is rejected rather than silently accepted or crashing with
# base R's unrelated "missing value where TRUE/FALSE needed".

test_that("sample() rejects a NaN draw from the density estimator instead of returning it", {
  set.seed(30)
  prior <- prior_uniform(0, 1)
  simulator <- function(theta) theta + stats::rnorm(1, sd = 0.05)
  fit <- npe(prior, simulator, n_simulations = 500,
             density_estimator = "linear_gaussian")
  post <- posterior(fit, x_obs = 0.5)

  # Corrupt only the first round's draw: a real density estimator would not
  # reliably reproduce a NaN row, so this pins down exactly one bad row and
  # lets later rounds recover normally, the same way rejecting an
  # out-of-bounds draw does.
  real_de_sample <- de_sample
  call_count <- 0L
  local_mocked_bindings(
    de_sample = function(de, x, n) {
      call_count <<- call_count + 1L
      draw <- real_de_sample(de, x, n)
      if (call_count == 1L) draw[1, ] <- NaN
      draw
    }
  )

  draws <- sample(post, n = 50, max_sampling_batches = 10)
  expect_equal(nrow(draws), 50L)
  expect_false(anyNA(draws))
  # without the fix, the NaN row is kept (not dropped), so round one alone
  # would already satisfy n and no second round would run
  expect_gt(call_count, 1L)
})

test_that("log_prob(normalize = TRUE) does not crash when normalization draws contain NaN", {
  set.seed(31)
  prior <- prior_uniform(0, 1)
  simulator <- function(theta) theta + stats::rnorm(1, sd = 0.05)
  fit <- npe(prior, simulator, n_simulations = 500,
             density_estimator = "linear_gaussian")
  post <- posterior(fit, x_obs = 0.5)

  real_de_sample <- de_sample
  local_mocked_bindings(
    de_sample = function(de, x, n) {
      draw <- real_de_sample(de, x, n)
      draw[1, ] <- NaN
      draw
    }
  )

  expect_no_error(lp <- log_prob(post, 0.5, n_normalization = 200))
  expect_true(is.finite(lp))
})

test_that("map_estimate()'s objective does not crash when queried at a non-finite point", {
  set.seed(32)
  prior <- prior_uniform(0, 1)
  simulator <- function(theta) theta + stats::rnorm(1, sd = 0.05)
  fit <- npe(prior, simulator, n_simulations = 500,
             density_estimator = "linear_gaussian")
  post <- posterior(fit, x_obs = 0.5)

  # Stand in for an optimizer step that proposes a non-finite parameter: call
  # the real objective at NaN before delegating to the real optim() so
  # map_estimate() still returns normally.
  queried <- NULL
  real_optim <- stats::optim
  local_mocked_bindings(
    optim = function(par, fn, ...) {
      queried <<- fn(NaN)
      real_optim(par, fn, ...)
    },
    .package = "stats"
  )

  expect_no_error(map <- map_estimate(post))
  expect_equal(queried, Inf)
  expect_true(within_support(prior, matrix(map, nrow = 1)))
})
