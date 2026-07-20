# Shared training engine: restarts, LR decay, gradient clipping, history.

test_that("training engine supports restarts and records history", {
  skip_if_no_torch()
  set.seed(3)
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + matrix(rnorm(length(theta), sd = 0.3),
                                              nrow = nrow(theta))
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
  simulator <- function(theta) theta + matrix(rnorm(length(theta), sd = 0.3),
                                              nrow = nrow(theta))
  fit <- npe(prior, simulator, n_simulations = 300, density_estimator = "mdn",
             n_components = 1L, hidden = c(10L), max_epochs = 10L,
             clip_grad_norm = Inf, seed = 4)
  expect_true(is.finite(fit$de$best_val_loss))
})
