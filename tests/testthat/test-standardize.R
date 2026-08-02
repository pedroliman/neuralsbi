# A constant column is left at scale 1, which is the only safe thing to do and
# also the reason it goes unnoticed: nothing downstream fails.

test_that("fit_standardizer() is silent unless asked to speak", {
  x <- cbind(a = rnorm(20), b = rep(3, 20))
  expect_silent(fit_standardizer(x))
  expect_equal(unname(fit_standardizer(x)$scale[2]), 1)
})

test_that("fit_standardizer() names constant columns and the side they are on", {
  x <- cbind(a = rnorm(20), flat = rep(3, 20))
  expect_warning(fit_standardizer(x, what = "x"),
                 "`x` has 1 constant column \\(flat\\)")
  expect_warning(fit_standardizer(x, what = "x"), "no information about `theta`")
  expect_warning(fit_standardizer(x, what = "theta"),
                 "`theta` has 1 constant column \\(flat\\)")
  expect_warning(fit_standardizer(x, what = "theta"), "cannot be identified")
})

test_that("fit_standardizer() falls back to the column index", {
  x <- cbind(rnorm(20), rep(3, 20), rep(-1, 20))
  expect_warning(fit_standardizer(x, what = "x"),
                 "2 constant columns \\(column 2, column 3\\)")
  # Half-named input names what it can.
  y <- cbind(rnorm(20), rep(3, 20))
  colnames(y) <- c("first", "")
  expect_warning(fit_standardizer(y, what = "x"), "\\(column 2\\)")
})

test_that("fit_standardizer() reports a single row as undefined, not constant", {
  expect_warning(fit_standardizer(matrix(0, 1, 2), what = "x"),
                 "`x` has one row, so the standard deviation of 2 columns")
})

test_that("a constant column still standardizes to zero", {
  x <- cbind(a = 1:20 + 0, b = rep(3, 20))
  std <- fit_standardizer(x)
  z <- apply_standardizer(std, x)
  expect_true(all(z[, 2] == 0))
  expect_equal(invert_standardizer(std, z), x)
})

test_that("npe() warns about a constant summary statistic", {
  prior <- prior_uniform(c(mu = -2), c(mu = 2))
  sim <- function(mu) c(y = mu + stats::rnorm(1, sd = 0.1), const = 1)
  expect_warning(npe(prior, sim, n_simulations = 50,
                     density_estimator = "linear_gaussian"),
                 "`x` has 1 constant column \\(const\\)")
})

test_that("standardize = FALSE builds its degenerate standardizer quietly", {
  prior <- prior_uniform(c(mu = -2), c(mu = 2))
  sim <- function(mu) c(y = mu + stats::rnorm(1, sd = 0.1))
  expect_silent(npe(prior, sim, n_simulations = 50,
                    density_estimator = "linear_gaussian",
                    standardize = FALSE))
})
