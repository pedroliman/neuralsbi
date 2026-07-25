# summary()/as.data.frame() ergonomics and plot_coverage(), torch-free via
# the linear_gaussian estimator.

fit_lg <- function() {
  set.seed(7)
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + matrix(rnorm(length(theta), sd = 0.5),
                                              nrow = nrow(theta))
  npe(prior, simulator, n_simulations = 1500,
      density_estimator = "linear_gaussian", seed = 7)
}

test_that("samples convert to data frames and summarize per parameter", {
  fit <- fit_lg()
  post <- posterior(fit, x_obs = c(1, -0.5))
  draws <- sample(post, 2000)

  df <- as.data.frame(draws)
  expect_s3_class(df, "data.frame")
  expect_equal(dim(df), c(2000L, 2L))
  expect_equal(names(df), c("theta1", "theta2"))

  s <- summary(draws)
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), 2L)
  expect_true(all(c("parameter", "mean", "sd", "q50") %in% names(s)))
  expect_true(all(s$q2.5 < s$q50 & s$q50 < s$q97.5))

  s2 <- summary(post, n = 500)
  expect_equal(nrow(s2), 2L)
})

test_that("summary of a fit returns training info invisibly", {
  fit <- fit_lg()
  info <- withVisible(summary(fit))
  expect_false(info$visible)
  expect_equal(info$value$density_estimator, "linear_gaussian")
  expect_equal(info$value$n_simulations, 1500L)
})

test_that("plot_coverage runs on an sbc result and returns coverage", {
  skip_if_no_ggplot2()
  fit <- fit_lg()
  simulator <- function(theta) theta + matrix(rnorm(length(theta), sd = 0.5),
                                              nrow = nrow(theta))
  res <- sbc(fit, simulator, n_sbc = 40L, n_posterior_samples = 200L, seed = 8)
  pdf(NULL)
  on.exit(dev.off())
  cov <- plot_coverage(res)
  expect_s3_class(cov, "data.frame")
  expect_equal(names(cov)[1], "nominal")
  # rough calibration for the exact estimator
  expect_lt(max(abs(cov$param1 - cov$nominal)), 0.3)
})
