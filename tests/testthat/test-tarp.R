# TARP (Lemos et al. 2023) is a joint calibration test. With the exact
# linear_gaussian estimator the ECP curve must sit on the diagonal; running the
# same fit against a simulator with inflated noise must pull it off.

test_that("tarp is calibrated for the exact linear-Gaussian posterior", {
  set.seed(7)
  sigma <- 0.5
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) {
    theta + matrix(rnorm(length(theta), sd = sigma), nrow = nrow(theta))
  }
  fit <- npe(prior, simulator, n_simulations = 3000,
             density_estimator = "linear_gaussian")

  res <- tarp(fit, simulator, n_tarp = 150L, n_posterior_samples = 300L,
              seed = 1)
  expect_s3_class(res, "nsbi_tarp")
  expect_length(res$coverage_values, 150L)
  expect_true(all(res$coverage_values >= 0 & res$coverage_values <= 1))
  # ECP within Monte-Carlo noise of the diagonal
  expect_lt(max(abs(res$ecp - res$levels)), 0.12)
  # coverage values themselves look uniform
  # ties from the 1/300 quantization only make the KS test conservative here
  p <- suppressWarnings(stats::ks.test(res$coverage_values, "punif")$p.value)
  expect_gt(p, 0.01)
})

test_that("tarp detects an overconfident posterior", {
  set.seed(8)
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  narrow_sim <- function(theta) {
    theta + matrix(rnorm(length(theta), sd = 0.3), nrow = nrow(theta))
  }
  wide_sim <- function(theta) {
    theta + matrix(rnorm(length(theta), sd = 1.5), nrow = nrow(theta))
  }
  # trained for sd 0.3, evaluated on data with sd 1.5: posterior far too narrow
  fit <- npe(prior, narrow_sim, n_simulations = 3000,
             density_estimator = "linear_gaussian")
  res <- tarp(fit, wide_sim, n_tarp = 150L, n_posterior_samples = 300L,
              seed = 2)
  expect_gt(max(abs(res$ecp - res$levels)), 0.15)
})

test_that("tarp supports prior reference points and prints", {
  set.seed(9)
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) {
    theta + matrix(rnorm(length(theta), sd = 0.5), nrow = nrow(theta))
  }
  fit <- npe(prior, simulator, n_simulations = 1000,
             density_estimator = "linear_gaussian")
  res <- tarp(fit, simulator, n_tarp = 50L, n_posterior_samples = 100L,
              references = "prior", seed = 3)
  expect_s3_class(res, "nsbi_tarp")
  expect_output(print(res), "nsbi_tarp")
})

test_that("plot_tarp runs and returns the ECP curve", {
  set.seed(10)
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) {
    theta + matrix(rnorm(length(theta), sd = 0.5), nrow = nrow(theta))
  }
  fit <- npe(prior, simulator, n_simulations = 1000,
             density_estimator = "linear_gaussian")
  res <- tarp(fit, simulator, n_tarp = 40L, n_posterior_samples = 100L,
              seed = 4)
  path <- tempfile(fileext = ".png")
  grDevices::png(path)
  out <- plot_tarp(res)
  grDevices::dev.off()
  expect_true(file.exists(path))
  expect_named(out, c("nominal", "ecp"))
  unlink(path)
})
