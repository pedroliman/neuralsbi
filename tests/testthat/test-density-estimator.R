# linear_gaussian is the oracle the rest of the suite is checked against, so
# its own contract is worth pinning down: the closed-form fit, and the shape
# checks its two methods share with the neural estimators.

test_that("fit_linear_gaussian() recovers the analytic conditional", {
  set.seed(11)
  n <- 4000; p <- 2; s <- 0.5
  theta <- matrix(rnorm(n * p), ncol = p)
  x <- theta + matrix(rnorm(n * p, sd = s), ncol = p)

  de <- fit_linear_gaussian(theta, x)

  # theta ~ N(0, I), x | theta ~ N(theta, s^2 I): the conditional mean is
  # x / (1 + s^2) and the conditional covariance is s^2 / (1 + s^2) * I.
  shrink <- 1 / (1 + s^2)
  expect_equal(de$B[1, ], rep(0, p), tolerance = 0.05)
  expect_equal(de$B[-1, ], diag(shrink, p), tolerance = 0.05)
  expect_equal(de$Sigma, diag(s^2 * shrink, p), tolerance = 0.05)

  x_obs <- c(1, -0.5)
  expect_equal(as.numeric(lingauss_mean(de, x_obs)), x_obs * shrink,
               tolerance = 0.05)
})

test_that("fit_linear_gaussian() recovers the coefficients of a noiseless model", {
  # With theta an exact linear function of x there is one right answer for B,
  # and the fit has to land on it rather than near it.
  set.seed(20)
  x <- matrix(stats::rnorm(600), ncol = 3)
  B_true <- rbind(c(0.5, -1), c(1, 0), c(-0.25, 2), c(3, 0.1))
  theta <- cbind(1, x) %*% B_true

  de <- fit_linear_gaussian(theta, x, ridge = 1e-10)

  expect_lt(max(abs(unname(de$B) - B_true)), 1e-6)
  # Nothing is left over, so the ridge is all that keeps Sigma positive
  # definite. It is measured against each target column's own variance, since
  # the residuals have no scale left to measure against.
  expect_equal(diag(de$Sigma), 1e-10 * apply(theta, 2, stats::var),
               tolerance = 1e-4)
  expect_true(all(is.finite(de$chol)))
})

test_that("the ridge keeps a rank-deficient design solvable", {
  set.seed(21)
  x <- matrix(stats::rnorm(300), ncol = 3)
  x[, 3] <- x[, 1]                     # exactly collinear, so X'X is singular
  theta <- x[, 1, drop = FALSE] + stats::rnorm(100, sd = 0.1)

  de <- fit_linear_gaussian(theta, x)
  expect_s3_class(de, "nsbi_de_lingauss")
  expect_true(all(is.finite(de$B)))
  # The default ridge is what saves it: without one, solve() gives up.
  expect_error(fit_linear_gaussian(theta, x, ridge = 0), "singular")
})

test_that("the residual covariance divides by the residual degrees of freedom", {
  set.seed(22)
  n <- 30
  x <- matrix(stats::rnorm(n * 2), ncol = 2)
  theta <- cbind(x[, 1] - 0.5 * x[, 2] + stats::rnorm(n, sd = 0.4))
  ridge <- 1e-6

  de <- fit_linear_gaussian(theta, x, ridge = ridge)

  X <- cbind(1, x)
  resid <- theta - X %*% de$B
  # The ridge is relative to the target's variance, not an absolute constant.
  bump <- ridge * apply(theta, 2, stats::var)
  expect_equal(de$Sigma, crossprod(resid) / (n - ncol(X)) + bump)
  # n is the maximum-likelihood denominator, and at n = 30 the two differ.
  expect_false(isTRUE(all.equal(de$Sigma, crossprod(resid) / n + bump)))
})

test_that("the ridge does not depend on the scale of the data", {
  # An absolute ridge is only negligible when the data happen to be O(1). This
  # target has variance 2.5e-07, four orders below the 1e-06 default, and used
  # to come back five times too wide -- which is what nle(standardize = FALSE)
  # on a simulator reporting small numbers was getting.
  set.seed(24)
  n <- 4000
  x <- matrix(stats::rnorm(n * 2), ncol = 2)
  clean <- cbind(x[, 1] - 0.5 * x[, 2])

  for (scale in c(1e-3, 1, 1e3)) {
    theta <- (clean + stats::rnorm(n, sd = 0.5)) * scale
    de <- fit_linear_gaussian(theta, x)
    # The truth is 0.25 * scale^2; the tolerance is the sampling error in a
    # variance from n draws, which is what is left once the ridge is scale-free.
    expect_equal(as.numeric(de$Sigma) / scale^2, 0.25, tolerance = 0.05,
                 label = paste("scale", scale))
    expect_equal(as.numeric(de$B[2, 1]) / scale, 1, tolerance = 0.05,
                 label = paste("scale", scale))
  }
})

test_that("the denominator floors at 1 when the design is wider than it is tall", {
  # n - ncol(X) is negative here. A negative denominator would flip the sign of
  # Sigma and take chol() down with it.
  set.seed(23)
  de <- fit_linear_gaussian(matrix(stats::rnorm(4), ncol = 2),
                            matrix(stats::rnorm(8), ncol = 4), ridge = 0.1)

  expect_true(all(is.finite(de$Sigma)))
  expect_true(all(diag(de$Sigma) > 0))
  expect_true(all(is.finite(de$chol)))
})

test_that("fit_linear_gaussian() records both dimensions", {
  de <- fit_linear_gaussian(matrix(rnorm(300), ncol = 3),
                            matrix(rnorm(500), ncol = 5))
  expect_equal(de$dim_theta, 3L)
  expect_equal(de$dim_x, 5L)
})

test_that("the linear-Gaussian methods reject an x of the wrong width", {
  set.seed(12)
  de <- fit_linear_gaussian(matrix(rnorm(200), ncol = 2),
                            matrix(rnorm(300), ncol = 3))
  theta <- matrix(0, nrow = 1, ncol = 2)

  # Before dim_x was recorded these reached the matrix product in
  # lingauss_mean() and came back as "non-conformable arguments".
  expect_error(de_log_prob(de, theta, matrix(0, nrow = 1, ncol = 4)),
               "Expected 3 columns but got 4")
  expect_error(de_sample(de, matrix(0, nrow = 1, ncol = 2), 5),
               "Expected 3 columns but got 2")
  expect_error(lingauss_mean(de, matrix(0, nrow = 2, ncol = 5)),
               "Expected 3 columns but got 5")
})

test_that("a bare x of length dim_x is one observation, not a column", {
  set.seed(13)
  de <- fit_linear_gaussian(matrix(rnorm(200), ncol = 2),
                            matrix(rnorm(300), ncol = 3))

  mu <- lingauss_mean(de, c(0.1, 0.2, 0.3))
  expect_equal(dim(mu), c(1L, 2L))
  expect_equal(mu, lingauss_mean(de, matrix(c(0.1, 0.2, 0.3), nrow = 1)))

  expect_length(de_log_prob(de, c(0, 0), c(0.1, 0.2, 0.3)), 1L)
  expect_equal(dim(de_sample(de, c(0.1, 0.2, 0.3), 4)), c(4L, 2L))
})

test_that("de_log_prob() is dmvnorm_chol() around the fitted conditional mean", {
  set.seed(24)
  theta <- matrix(stats::rnorm(600), ncol = 2)
  x <- theta + matrix(stats::rnorm(600, sd = 0.5), ncol = 2)
  de <- fit_linear_gaussian(theta, x)

  x_obs <- matrix(c(0.7, -1.1), nrow = 1)
  mu <- as.numeric(lingauss_mean(de, x_obs))
  query <- rbind(mu, mu + c(0.3, -0.2), mu - c(1, 1))

  lp <- de_log_prob(de, query, x_obs)
  expect_equal(lp, dmvnorm_chol(query, mu, de$chol))
  # The conditional mean is the mode, so the first row scores highest.
  expect_equal(which.max(lp), 1L)
})

test_that("a single x is broadcast across many theta rows", {
  set.seed(25)
  theta <- matrix(stats::rnorm(400), ncol = 2)
  x <- theta + matrix(stats::rnorm(400, sd = 0.5), ncol = 2)
  de <- fit_linear_gaussian(theta, x)

  query <- matrix(stats::rnorm(20), ncol = 2)
  x_obs <- matrix(c(0.4, -0.9), nrow = 1)
  repeated <- matrix(c(0.4, -0.9), nrow = nrow(query), ncol = 2, byrow = TRUE)

  expect_length(de_log_prob(de, query, x_obs), 10L)
  expect_equal(de_log_prob(de, query, x_obs), de_log_prob(de, query, repeated))
})

test_that("de_sample.nsbi_de_lingauss() draws from the fitted conditional", {
  set.seed(14)
  n <- 4000; s <- 0.5
  theta <- matrix(rnorm(n * 2), ncol = 2)
  x <- theta + matrix(rnorm(n * 2, sd = s), ncol = 2)
  de <- fit_linear_gaussian(theta, x)

  x_obs <- c(1, -0.5)
  draws <- de_sample(de, x_obs, 20000)
  expect_equal(colMeans(draws), as.numeric(lingauss_mean(de, x_obs)),
               tolerance = 0.05)
  expect_equal(stats::cov(draws), de$Sigma, tolerance = 0.05)
})

test_that("an estimator fitted before dim_x existed still evaluates", {
  # Serialized fits carry no dim_x, and de_rebuild_net() has no lingauss
  # branch to add one, so the methods must not require it.
  set.seed(15)
  de <- fit_linear_gaussian(matrix(rnorm(200), ncol = 2),
                            matrix(rnorm(300), ncol = 3))
  old <- de
  old$dim_x <- NULL
  expect_equal(lingauss_mean(old, matrix(c(0.1, 0.2, 0.3), nrow = 1)),
               lingauss_mean(de, matrix(c(0.1, 0.2, 0.3), nrow = 1)))
})
