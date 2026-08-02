# Argument guards on the plotting functions. Drawing needs ggplot2, and
# pairplot() also needs GGally and ggdensity; all three are Suggests, so these
# use the skips from helper-suggests.R.

sbc_fixture <- function(seed = 1) {
  set.seed(seed)
  prior <- prior_normal(mean = c(mu = 0, nu = 0), sd = 1)
  simulator <- function(mu, nu) c(a = mu, b = nu) + stats::rnorm(2, sd = 0.4)
  fit <- npe(prior, simulator, n_simulations = 600,
             density_estimator = "linear_gaussian")
  # 100 trials over the 20 rank bins keeps chisq.test() off its
  # small-expected-count warning.
  sbc(fit, simulator, n_sbc = 100L, n_posterior_samples = 50L, seed = seed)
}

test_that("plot_sbc() refuses a parameter it cannot subscript", {
  skip_if_no_ggplot2()
  res <- sbc_fixture()

  expect_error(plot_sbc(res, param = 99),
               regexp = "`param` must be one parameter index between 1 and 2")
  # A name is an accepted value, so the names belong in the message.
  expect_error(plot_sbc(res, param = 99), regexp = "or one of mu, nu")
  expect_error(plot_sbc(res, param = 0), regexp = "between 1 and 2")
  expect_error(plot_sbc(res, param = 1.5), regexp = "between 1 and 2")
  expect_error(plot_sbc(res, param = c(1, 2)),
               regexp = "not a length-2 numeric vector")
  expect_error(plot_sbc(res, param = NA), regexp = "not NA")
})

test_that("plot_sbc() resolves a parameter name against the rank columns", {
  skip_if_no_ggplot2()
  res <- sbc_fixture()

  expect_error(plot_sbc(res, param = "sigma"),
               regexp = "not one of the parameter names: mu, nu")

  unnamed <- res
  colnames(unnamed$ranks) <- NULL
  expect_error(plot_sbc(unnamed, param = "mu"),
               regexp = "the parameters are unnamed")

  path <- tempfile(fileext = ".png")
  grDevices::png(path)
  by_name <- plot_sbc(res, param = "nu")
  by_index <- plot_sbc(res, param = 2L)
  grDevices::dev.off()

  expect_identical(by_name, by_index)
  expect_length(by_name, 100L)
  expect_true(file.exists(path))
  unlink(path)
})

test_that("pairplot() checks limits against the number of parameters", {
  skip_if_no_ggally()
  set.seed(2)
  draws <- matrix(stats::rnorm(1500), ncol = 3)

  expect_error(pairplot(draws, limits = list(c(-3, 3), c(-3, 3))),
               regexp = "`limits` has 2 elements but `samples` has 3 parameters")
  expect_error(pairplot(draws, limits = matrix(c(-3, 3), nrow = 1)),
               regexp = "`limits` has 1 row but `samples` has 3 parameters")
  expect_error(pairplot(draws, limits = c(-3, 3)),
               regexp = "`limits` must be a list of 3 elements")
})

test_that("pairplot() accepts one limit pair per parameter", {
  skip_if_no_ggally()
  set.seed(3)
  draws <- matrix(stats::rnorm(1500), ncol = 3)
  path <- tempfile(fileext = ".png")

  grDevices::png(path)
  from_list <- pairplot(draws, limits = list(c(-3, 3), c(-3, 3), c(-3, 3)))
  from_matrix <- pairplot(draws, limits = matrix(rep(c(-3, 3), each = 3),
                                                 nrow = 3))
  grDevices::dev.off()

  expect_s3_class(from_list, "ggmatrix")
  expect_s3_class(from_matrix, "ggmatrix")
  expect_true(file.exists(path))
  unlink(path)
})
