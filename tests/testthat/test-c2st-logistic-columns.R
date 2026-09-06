# c2st_logistic_prob() built the training frame as data.frame(y = y_train,
# x_train) and the prediction frame as data.frame(x_test). data.frame() names
# a lone unnamed column after the deparsed argument it was given, so a
# single-column, unnamed x_train/x_test pair got mismatched column names
# ("x_train" vs "x_test"), glm() fit a coefficient under one name, and
# predict() errored looking for it under the other (#278).

test_that("c2st() with classifier = \"logistic\" handles single-column, unnamed matrices", {
  set.seed(1)
  a <- matrix(rnorm(300), ncol = 1)
  b <- matrix(rnorm(300, 3), ncol = 1)
  expect_null(colnames(a))
  expect_null(colnames(b))

  out <- c2st(a, b, classifier = "logistic", seed = 1)
  expect_true(is.finite(out$accuracy))
  expect_true(is.finite(out$auc))
  # a and b are well separated (mean shift of 3), so a linear classifier
  # should tell them apart easily.
  expect_gt(out$accuracy, 0.8)
})

test_that("c2st() with classifier = \"logistic\" still handles multi-column matrices", {
  set.seed(2)
  a <- matrix(rnorm(600), ncol = 2)
  b <- cbind(rnorm(300, mean = 3), rnorm(300, mean = -3))
  expect_null(colnames(a))
  expect_null(colnames(b))

  out <- c2st(a, b, classifier = "logistic", seed = 1)
  expect_true(is.finite(out$accuracy))
  expect_gt(out$accuracy, 0.8)
})

test_that("c2st() with classifier = \"logistic\" ignores pre-existing column names", {
  set.seed(3)
  a <- matrix(rnorm(300), ncol = 1, dimnames = list(NULL, "theta"))
  b <- matrix(rnorm(300, 3), ncol = 1, dimnames = list(NULL, "other_name"))

  out <- c2st(a, b, classifier = "logistic", seed = 1)
  expect_true(is.finite(out$accuracy))
  expect_gt(out$accuracy, 0.8)
})
