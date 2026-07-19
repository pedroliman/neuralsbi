test_that("uniform prior samples within bounds and has correct shape", {
  prior <- prior_uniform(low = c(-2, 0), high = c(2, 5))
  th <- sample_prior(prior, 100)
  expect_equal(dim(th), c(100L, 2L))
  expect_true(all(th[, 1] >= -2 & th[, 1] <= 2))
  expect_true(all(th[, 2] >= 0 & th[, 2] <= 5))
})

test_that("uniform prior log_prob is constant inside, -Inf outside", {
  prior <- prior_uniform(low = c(0, 0), high = c(1, 1))
  inside <- matrix(c(0.5, 0.5), nrow = 1)
  outside <- matrix(c(1.5, 0.5), nrow = 1)
  expect_equal(prior$log_prob(inside), 0)          # log(1/(1*1)) = 0
  expect_equal(prior$log_prob(outside), -Inf)
})

test_that("normal prior log_prob matches dnorm", {
  prior <- prior_normal(mean = c(0, 1), sd = c(1, 2))
  th <- matrix(c(0, 1), nrow = 1)
  expect_equal(prior$log_prob(th),
               dnorm(0, 0, 1, log = TRUE) + dnorm(1, 1, 2, log = TRUE))
})

test_that("within_support flags out-of-bounds rows", {
  prior <- prior_uniform(c(-1, -1), c(1, 1))
  th <- rbind(c(0, 0), c(2, 0), c(0, -3))
  expect_equal(within_support(prior, th), c(TRUE, FALSE, FALSE))
})

test_that("unbounded prior treats everything as in-support", {
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  expect_true(all(within_support(prior, matrix(rnorm(20), ncol = 2))))
})
