# Shared training engine: restarts, LR decay, gradient clipping, history.

test_that("training engine supports restarts and records history", {
  skip_if_no_torch()
  set.seed(3)
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + rnorm(length(theta), sd = 0.3)
  fit <- npe(prior, simulator, n_simulations = 500, density_estimator = "mdn",
             n_components = 1L, hidden = c(20L), max_epochs = 30L,
             n_restarts = 2L, seed = 3)
  expect_true(is.finite(fit$de$best_val_loss))
  hist <- fit$de$history
  expect_s3_class(hist, "data.frame")
  expect_named(hist, c("epoch", "train_loss", "val_loss"))
  expect_gt(nrow(hist), 0)
  expect_true(all(is.finite(hist$val_loss)))
})

test_that("clip_grad_norm = Inf disables clipping without breaking training", {
  skip_if_no_torch()
  set.seed(4)
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + rnorm(length(theta), sd = 0.3)
  fit <- npe(prior, simulator, n_simulations = 300, density_estimator = "mdn",
             n_components = 1L, hidden = c(10L), max_epochs = 10L,
             clip_grad_norm = Inf, seed = 4)
  expect_true(is.finite(fit$de$best_val_loss))
})

test_that("train_conditional_de() checks its controls without torch", {
  theta <- matrix(stats::rnorm(50), ncol = 1)
  x <- matrix(stats::rnorm(50), ncol = 1)
  # build_net and log_prob_fn are never reached: the controls are checked
  # before the loop, which is the point of checking them there.
  train <- function(...) {
    train_conditional_de(build_net = function() stop("not reached"),
                         log_prob_fn = function(...) stop("not reached"),
                         theta = theta, x = x, ...)
  }
  expect_error(train(batch_size = 0), "`batch_size` must be")
  expect_error(train(max_epochs = 0), "`max_epochs` must be")
  expect_error(train(patience = -1), "`patience` must be")
  expect_error(train(n_restarts = 0), "`n_restarts` must be")
  expect_error(train(lr = -1e-3), "`lr` must be a single positive number")
  expect_error(train(clip_grad_norm = 0),
               "`clip_grad_norm` must be a single positive number or Inf")
  expect_error(train(validation_fraction = 1), "strictly between 0 and 1")
  expect_error(train(validation_fraction = 0), "strictly between 0 and 1")
})

test_that("a torch-backed fit records the device it trained on", {
  skip_if_no_torch()
  set.seed(5)
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + rnorm(length(theta), sd = 0.3)
  fit <- npe(prior, simulator, n_simulations = 200, density_estimator = "mdn",
             n_components = 1L, hidden = c(8L), max_epochs = 5L, seed = 5)
  # No GPU on the CI runner or in this sandbox, so "cpu" is the only device
  # actually reachable here; the cuda/mps fallback path is exercised by
  # resolve_device()'s own unit tests instead.
  expect_identical(fit$de$device, "cpu")
})

test_that("train_conditional_de() rejects a malformed device before training", {
  theta <- matrix(stats::rnorm(50), ncol = 1)
  x <- matrix(stats::rnorm(50), ncol = 1)
  expect_error(
    train_conditional_de(build_net = function() stop("not reached"),
                         log_prob_fn = function(...) stop("not reached"),
                         theta = theta, x = x, device = "gpu"),
    "`device` must be one of")
})

test_that("train_conditional_de() says how many rows the split would need", {
  one <- matrix(0, nrow = 1, ncol = 1)
  expect_error(
    train_conditional_de(build_net = function() stop("not reached"),
                         log_prob_fn = function(...) stop("not reached"),
                         theta = one, x = one, validation_fraction = 0.5),
    "holds out 1 row of 1, leaving nothing to train on\\. At this fraction the estimator needs at least 2 rows")
})
