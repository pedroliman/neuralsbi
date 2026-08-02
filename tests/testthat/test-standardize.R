# Every estimator trains in standardized space, so the round trip is the one
# property the whole design rests on: whatever goes through apply_standardizer()
# has to come back through invert_standardizer() unchanged.

test_that("a round trip through the standardizer returns the input", {
  set.seed(101)
  x <- cbind(big = stats::rnorm(50, mean = 500, sd = 40),
             small = stats::rnorm(50, mean = -0.02, sd = 0.005),
             plain = stats::runif(50))
  std <- fit_standardizer(x)
  z <- apply_standardizer(std, x)

  # Absolute error, not relative: the centered values pass through zero, where
  # a relative comparison says nothing useful.
  expect_lt(max(abs(invert_standardizer(std, z) - x)), 1e-10)
  expect_lt(max(abs(colMeans(z))), 1e-10)
  expect_lt(max(abs(apply(z, 2, stats::sd) - 1)), 1e-10)
})

test_that("invert_standardizer() holds for a z the standardizer never produced", {
  # Draws arrive from the estimator in standardized space rather than from
  # apply_standardizer(), so the inverse has to hold for an arbitrary z.
  std <- fit_standardizer(cbind(a = c(1, 3, 5, 7), b = c(-2, 0, 2, 4)))
  z <- cbind(c(-1.5, 0, 2.25), c(0.5, -3, 1))

  expect_lt(max(abs(apply_standardizer(std, invert_standardizer(std, z)) - z)),
            1e-10)
})

test_that("standardizer_log_jac() is minus the summed log scale", {
  x <- cbind(a = c(1, 2, 3, 4, 10), b = c(0, 5, -5, 2, 1))
  std <- fit_standardizer(x)

  expect_equal(standardizer_log_jac(std), -sum(log(apply(x, 2, stats::sd))))
  # A standardizer that scales nothing contributes nothing to the density.
  expect_equal(standardizer_log_jac(fit_standardizer(matrix(0, 1, 3))), 0)
})

test_that("the no-spread guard leaves scale at 1 in both branches", {
  # sd() is 0 for a constant column and NA for a single row. Dividing by either
  # gives Inf or NaN, so both fall back to scale 1.
  flat <- fit_standardizer(cbind(a = stats::rnorm(20), b = rep(3, 20)))
  expect_equal(unname(flat$scale[2]), 1)
  expect_true(all(is.finite(flat$scale)))

  one_row <- fit_standardizer(matrix(c(4, -7), nrow = 1))
  expect_equal(unname(one_row$scale), c(1, 1))
  expect_equal(unname(one_row$center), c(4, -7))
  expect_equal(unname(apply_standardizer(one_row, matrix(c(4, -7), nrow = 1))),
               matrix(0, 1, 2))
})

test_that("eps decides how little spread counts as none", {
  tiny <- cbind(a = rep(c(0, 1e-12), 10))
  expect_equal(unname(fit_standardizer(tiny)$scale), 1)
  expect_lt(unname(fit_standardizer(tiny, eps = 1e-20)$scale), 1e-10)
})

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
