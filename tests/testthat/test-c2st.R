# z_score was only ever tested with isTRUE(), so "TRUE" or 1 silently took
# the "don't z-score" branch instead of erroring, training the classifier on
# the raw scale of x and y with nothing said about it (#327).

test_that("c2st() rejects a non-logical z_score instead of silently skipping it", {
  a <- matrix(rnorm(400), ncol = 2)
  b <- matrix(rnorm(400), ncol = 2)

  expect_error(c2st(a, b, classifier = "logistic", seed = 1, z_score = 1),
               "`z_score` must be TRUE or FALSE")
  expect_error(c2st(a, b, classifier = "logistic", seed = 1, z_score = "TRUE"),
               "`z_score` must be TRUE or FALSE")
})

test_that("c2st() still runs with z_score = TRUE and z_score = FALSE", {
  a <- matrix(rnorm(400), ncol = 2)
  b <- matrix(rnorm(400), ncol = 2)

  out_true <- c2st(a, b, classifier = "logistic", seed = 1, z_score = TRUE)
  expect_true(is.finite(out_true$accuracy))

  out_false <- c2st(a, b, classifier = "logistic", seed = 1, z_score = FALSE)
  expect_true(is.finite(out_false$accuracy))
})
