# A single NA/NaN/Inf anywhere in x or y used to standardize into a whole
# corrupted column (fit_standardizer()/apply_standardizer() carry no na.rm),
# and land several frames away in glm() or the MLP training loop with an error
# that names neither x nor y. c2st() should reject the entry itself, before
# any of that runs.

test_that("c2st() rejects a non-finite entry in x and names it", {
  a <- matrix(rnorm(20), ncol = 2)
  b <- matrix(rnorm(20), ncol = 2)

  a_na <- a
  a_na[3, 1] <- NA
  expect_error(c2st(a_na, b, classifier = "logistic", seed = 1),
               "`x` contains 1 non-finite value \\(NA\\), first at row 3, column 1")

  a_nan <- a
  a_nan[4, 2] <- NaN
  expect_error(c2st(a_nan, b, classifier = "logistic", seed = 1),
               "`x` contains 1 non-finite value \\(NaN\\)")

  a_inf <- a
  a_inf[1, 1] <- Inf
  expect_error(c2st(a_inf, b, classifier = "logistic", seed = 1),
               "`x` contains 1 non-finite value \\(Inf\\)")
})

test_that("c2st() rejects a non-finite entry in y and names it", {
  a <- matrix(rnorm(20), ncol = 2)
  b <- matrix(rnorm(20), ncol = 2)

  b_na <- b
  b_na[5, 2] <- NA
  expect_error(c2st(a, b_na, classifier = "logistic", seed = 1),
               "`y` contains 1 non-finite value \\(NA\\), first at row 5, column 2")
})

test_that("c2st() still runs on fully finite draws", {
  a <- matrix(rnorm(400), ncol = 2)
  b <- matrix(rnorm(400), ncol = 2)
  out <- c2st(a, b, classifier = "logistic", seed = 1)
  expect_true(is.finite(out$accuracy))
})
