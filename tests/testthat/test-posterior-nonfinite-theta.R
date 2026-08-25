# A NaN/NA entry in theta used to slip past log_prob.nsbi_posterior()'s
# validation: within_support() on a NaN row returns logical NA, and
# `lp[!within_support(prior, theta)] <- -Inf` is a no-op for an NA index, so
# the row fell through to de_log_prob() and came back as a silent NaN
# log-density instead of a named error (#221).

npe_fit_2d <- function(seed = 1) {
  set.seed(seed)
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + stats::rnorm(2, sd = 0.5)
  npe(prior, simulator, n_simulations = 400,
      density_estimator = "linear_gaussian")
}

test_that("log_prob() rejects a non-finite theta and names the entry", {
  fit <- npe_fit_2d()
  post <- posterior(fit, x_obs = c(0.4, -0.4))

  expect_error(log_prob(post, c(0.5, NA)),
               "`theta` contains 1 non-finite value \\(NA\\), first at position 2")
  expect_error(log_prob(post, c(NaN, 0.5)), "`theta` contains")
  expect_error(log_prob(post, matrix(c(1, 2, 3, NA), ncol = 2)),
               "`theta` contains 1 non-finite value \\(NA\\), first at row 2, column 2")
})

test_that("log_prob() still allows an Inf theta", {
  fit <- npe_fit_2d()
  post <- posterior(fit, x_obs = c(0.4, -0.4))

  # Inf is not an error: it resolves to -Inf log-density via the bounded
  # prior's out-of-support masking, or evaluates normally for an unbounded
  # one -- it is only NA/NaN that has no well-defined meaning.
  lp <- log_prob(post, c(0.5, Inf))
  expect_true(is.finite(lp) || lp == -Inf)
  expect_false(is.nan(lp))
})
