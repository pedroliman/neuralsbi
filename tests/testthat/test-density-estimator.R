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
