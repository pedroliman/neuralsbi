# Argument checking at the npe()/nle()/simulate_for_sbi() boundary. Every one
# of these calls used to reach base R with a bad number, usually after paying
# for the simulations first.

toy_prior <- function() prior_uniform(c(mu = -2), c(mu = 2))
toy_simulator <- function(mu) c(y = mu + stats::rnorm(1, sd = 0.1))

# A simulator that records every call, to prove the budget is not spent before
# the arguments are checked.
counting_simulator <- function(env) {
  function(mu) {
    env$calls <- env$calls + 1L
    c(y = mu + stats::rnorm(1, sd = 0.1))
  }
}

test_that("npe() rejects a simulation budget it cannot fit", {
  expect_error(npe(toy_prior(), toy_simulator, n_simulations = 0,
                   density_estimator = "linear_gaussian"),
               "`n_simulations` must be a single whole number of at least 2")
  expect_error(npe(toy_prior(), toy_simulator, n_simulations = -5,
                   density_estimator = "linear_gaussian"),
               "`n_simulations` must be")
  expect_error(npe(toy_prior(), toy_simulator, n_simulations = 1,
                   density_estimator = "linear_gaussian"),
               "not 1")
  expect_error(npe(toy_prior(), toy_simulator, n_simulations = 500.5,
                   density_estimator = "linear_gaussian"),
               "not 500.5")
})

test_that("npe() checks the controls before running the simulator", {
  env <- new.env(parent = emptyenv())
  env$calls <- 0L
  sim <- counting_simulator(env)

  expect_error(npe(toy_prior(), sim, n_simulations = 100,
                   density_estimator = "linear_gaussian",
                   validation_fraction = 1),
               "`validation_fraction` must be a single number strictly")
  expect_error(npe(toy_prior(), sim, n_simulations = 100,
                   density_estimator = "linear_gaussian", batch_size = 0),
               "`batch_size` must be")
  expect_error(npe(toy_prior(), sim, n_simulations = 100,
                   density_estimator = "linear_gaussian", n_restarts = 0),
               "`n_restarts` must be a single whole number of at least 1 since")
  expect_error(npe(toy_prior(), sim, n_simulations = 100,
                   density_estimator = "linear_gaussian", lr = 0),
               "`lr` must be a single positive number")
  expect_error(npe(toy_prior(), sim, n_simulations = 100,
                   density_estimator = "linear_gaussian", max_epochs = NA),
               "`max_epochs` must be")
  expect_error(npe(toy_prior(), sim, n_simulations = 100,
                   density_estimator = "linear_gaussian", patience = "20"),
               "`patience` must be")
  expect_error(npe(toy_prior(), sim, n_simulations = 100,
                   density_estimator = "linear_gaussian", clip_grad_norm = 0),
               "`clip_grad_norm` must be a single positive number or Inf")
  expect_identical(env$calls, 0L)
})

test_that("npe() checks architecture arguments", {
  expect_error(npe(toy_prior(), toy_simulator, n_simulations = 100,
                   density_estimator = "linear_gaussian", n_components = 0),
               "`n_components` must be")
  expect_error(npe(toy_prior(), toy_simulator, n_simulations = 100,
                   density_estimator = "linear_gaussian", n_transforms = 2.5),
               "`n_transforms` must be")
  expect_error(npe(toy_prior(), toy_simulator, n_simulations = 100,
                   density_estimator = "linear_gaussian", hidden = c(50, 0)),
               "`hidden` must be whole numbers of at least 1")
  expect_error(npe(toy_prior(), toy_simulator, n_simulations = 100,
                   density_estimator = "linear_gaussian", hidden = numeric(0)),
               "`hidden` must be")
  expect_error(npe(toy_prior(), toy_simulator, n_simulations = 100,
                   density_estimator = "linear_gaussian", n_bins = -1),
               "`n_bins` must be")
  expect_error(npe(toy_prior(), toy_simulator, n_simulations = 100,
                   density_estimator = "linear_gaussian", tail_bound = 0),
               "`tail_bound` must be a single positive number")
})

test_that("npe() names the prior when it is not one", {
  expect_error(npe(list(dim = 1), toy_simulator, n_simulations = 100,
                   density_estimator = "linear_gaussian"),
               "`prior` must be an nsbi_prior object, not list")
})

test_that("npe() rejects an unknown density estimator before simulating", {
  env <- new.env(parent = emptyenv())
  env$calls <- 0L
  expect_error(npe(toy_prior(), counting_simulator(env), n_simulations = 100,
                   density_estimator = "spline"),
               "'arg' should be one of")
  expect_identical(env$calls, 0L)
})

test_that("npe() leaves a valid call alone", {
  fit <- npe(toy_prior(), toy_simulator, n_simulations = 200,
             density_estimator = "linear_gaussian", seed = 1)
  expect_s3_class(fit, "nsbi_npe")
  expect_identical(fit$n_simulations, 200L)

  # n_simulations is ignored when the simulations are supplied, so a nonsense
  # value there is not the user's problem.
  theta <- matrix(stats::runif(100, -2, 2), ncol = 1)
  x <- theta + stats::rnorm(100, sd = 0.1)
  expect_s3_class(npe(toy_prior(), theta = theta, x = x, n_simulations = 0,
                      density_estimator = "linear_gaussian"),
                  "nsbi_npe")
})

test_that("npe() rejects simulations too few to split", {
  theta <- matrix(0.5, ncol = 1)
  x <- matrix(0.4, ncol = 1)
  expect_error(npe(toy_prior(), theta = theta, x = x,
                   density_estimator = "mdn"),
               "holds out 1 row of 1, leaving nothing to train on")
})

test_that("nle() checks the same controls", {
  expect_error(nle(toy_prior(), toy_simulator, n_simulations = 0,
                   density_estimator = "linear_gaussian"),
               "`n_simulations` must be")
  expect_error(nle(toy_prior(), toy_simulator, n_simulations = 100,
                   density_estimator = "linear_gaussian", batch_size = -1),
               "`batch_size` must be")
  expect_error(nle(toy_prior(), toy_simulator, n_simulations = 100,
                   density_estimator = "linear_gaussian", hidden = NULL),
               "`hidden` must be")
  expect_error(nle(list(), toy_simulator, n_simulations = 100,
                   density_estimator = "linear_gaussian"),
               "`prior` must be an nsbi_prior object")
})

test_that("simulate_for_sbi() checks n and the prior", {
  expect_error(simulate_for_sbi(toy_simulator, toy_prior(), n = 0),
               "`n` must be a single whole number of at least 1")
  expect_error(simulate_for_sbi(toy_simulator, toy_prior(), n = c(10, 20)),
               "a length-2 numeric vector")
  expect_error(simulate_for_sbi(toy_simulator, prior = NULL, n = 10),
               "`prior` must be an nsbi_prior object")
  expect_equal(nrow(simulate_for_sbi(toy_simulator, toy_prior(), n = 5)$theta),
               5L)
})
