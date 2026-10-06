# GitHub #381: a bare set.seed(seed) in the entry points left the caller's
# base R RNG in a state fixed by `seed`. Each test checks that the caller's
# next draw matches the draw from an untouched stream.

baseline_draw <- function() {
  set.seed(10)
  stats::runif(1)
}

seeded_fit <- function() {
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + stats::rnorm(length(theta), sd = 0.3)
  list(prior = prior, simulator = simulator,
       fit = npe(prior, simulator, n_simulations = 60,
                 density_estimator = "linear_gaussian", seed = 1))
}

test_that("local_seed() seeds the frame and restores the caller's stream", {
  f <- function() {
    local_seed(5)
    stats::runif(1)
  }
  set.seed(10)
  inside <- f()
  expect_identical(stats::runif(1), baseline_draw())
  set.seed(5)
  expect_identical(inside, stats::runif(1))
  expect_identical(f(), inside)
})

test_that("local_seed() removes .Random.seed when there was none", {
  f <- function() local_seed(5)
  if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
  f()
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("npe() with a seed leaves the caller's RNG stream alone", {
  set.seed(10)
  s <- seeded_fit()
  expect_identical(stats::runif(1), baseline_draw())
})

test_that("simulate_for_sbi() with a seed leaves the caller's RNG stream alone", {
  set.seed(10)
  simulate_for_sbi(function(theta) theta + stats::rnorm(length(theta)),
                   prior_normal(0, 1), n = 20, seed = 2)
  expect_identical(stats::runif(1), baseline_draw())
})

test_that("sbc() and tarp() with a seed leave the caller's RNG stream alone", {
  s <- seeded_fit()
  set.seed(10)
  sbc(s$fit, s$simulator, n_sbc = 20, n_posterior_samples = 50, seed = 3)
  expect_identical(stats::runif(1), baseline_draw())
  set.seed(10)
  tarp(s$fit, s$simulator, n_tarp = 20, n_posterior_samples = 50, seed = 3)
  expect_identical(stats::runif(1), baseline_draw())
})

test_that("c2st() with a seed leaves the caller's RNG stream alone", {
  a <- matrix(stats::rnorm(200), ncol = 2)
  b <- matrix(stats::rnorm(200), ncol = 2)
  set.seed(10)
  c2st(a, b, classifier = "logistic", seed = 4)
  expect_identical(stats::runif(1), baseline_draw())
})

test_that("npe_sequential() with a seed leaves the caller's RNG stream alone", {
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + stats::rnorm(length(theta), sd = 0.3)
  set.seed(10)
  npe_sequential(prior, simulator, x_obs = 0.2, n_rounds = 2L,
                 n_simulations = 60, density_estimator = "linear_gaussian",
                 seed = 5)
  expect_identical(stats::runif(1), baseline_draw())
})

test_that("nle() and nre() with a seed leave the caller's RNG stream alone", {
  skip_if_no_torch()
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + stats::rnorm(length(theta), sd = 0.3)
  set.seed(10)
  nle(prior, simulator, n_simulations = 100, density_estimator = "mdn",
      n_components = 1L, hidden = 8L, max_epochs = 3L, seed = 6)
  expect_identical(stats::runif(1), baseline_draw())
  set.seed(10)
  nre(prior, simulator, n_simulations = 100, hidden = 8L, max_epochs = 3L,
      seed = 6)
  expect_identical(stats::runif(1), baseline_draw())
})
