# GitHub #275: train_restarts() and c2st()'s MLP path both called
# torch::torch_manual_seed(seed) directly whenever a caller passed a seed, and
# never restored torch's prior global RNG state afterwards -- the torch
# analogue of #272's R-RNG leak in surrogate_potential(). Since torch has one
# global generator with no way to ask for an independent stream, a single
# seeded npe()/nle()/nre() fit or a seeded c2st() call left every later
# *unseeded* torch call in the same session drawing from wherever the seeded
# call left the generator, instead of from a fresh stream.

test_that("a seeded npe() fit does not mutate torch's global RNG stream", {
  skip_if_no_torch()
  torch::torch_manual_seed(123)
  before <- torch::as_array(torch::torch_get_rng_state())

  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + rnorm(length(theta), sd = 0.3)
  npe(prior, simulator, n_simulations = 200, density_estimator = "mdn",
      n_components = 1L, hidden = c(10L), max_epochs = 5L, seed = 99)

  after <- torch::as_array(torch::torch_get_rng_state())
  expect_identical(after, before)
})

test_that("a seeded c2st() MLP call does not mutate torch's global RNG stream", {
  skip_if_no_torch()
  torch::torch_manual_seed(321)
  before <- torch::as_array(torch::torch_get_rng_state())

  a <- matrix(stats::rnorm(200), ncol = 2)
  b <- matrix(stats::rnorm(200), ncol = 2)
  c2st(a, b, classifier = "mlp", seed = 5, max_epochs = 5L)

  after <- torch::as_array(torch::torch_get_rng_state())
  expect_identical(after, before)
})

test_that("set_torch_seed() returns the state needed to restore torch's RNG", {
  skip_if_no_torch()
  torch::torch_manual_seed(7)
  before <- torch::as_array(torch::torch_get_rng_state())

  old <- set_torch_seed(42)
  torch::torch_randn(3)  # spend some of the newly seeded stream
  expect_false(identical(torch::as_array(torch::torch_get_rng_state()), before))

  torch::torch_set_rng_state(old)
  expect_identical(torch::as_array(torch::torch_get_rng_state()), before)
})
