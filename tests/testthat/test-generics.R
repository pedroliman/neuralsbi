test_that("sample_posterior() dispatches on the posterior's class", {
  # sample_posterior() used to call sample.nsbi_posterior() rather than the
  # generic. An nle() posterior inherits nsbi_posterior, so the NPE forward-pass
  # sampler ran against an estimator with theta and x swapped: an error here,
  # where dim_theta is 2 and dim_x is 1, and silently wrong draws wherever the
  # two dimensions happen to match.
  set.seed(40)
  prior <- prior_uniform(c(mu = -3, nu = -3), c(mu = 3, nu = 3))
  simulator <- function(mu, nu) c(y = mu + 0.5 * nu + stats::rnorm(1, sd = 0.3))
  fit <- nle(prior, simulator, n_simulations = 2000,
             density_estimator = "linear_gaussian", seed = 41)
  x_obs <- matrix(stats::rnorm(30, mean = 1, sd = 0.3), ncol = 1)
  post <- posterior(fit, x_obs, n_chains = 4, warmup = 50, thin = 2, seed = 42)

  draws <- sample_posterior(post, n = 400)

  expect_equal(dim(draws), c(400L, 2L))
  expect_equal(colnames(draws), c("mu", "nu"))
  expect_equal(draws, sample(post, 400))
})

test_that("sample_posterior() agrees with sample() on an NPE posterior", {
  set.seed(43)
  prior <- prior_uniform(c(mu = -3), c(mu = 3))
  fit <- npe(prior, function(mu) c(y = stats::rnorm(1, mu, 0.5)),
             n_simulations = 1000, density_estimator = "linear_gaussian",
             seed = 44)
  post <- posterior(fit, matrix(1, nrow = 1))

  set.seed(45)
  direct <- sample(post, 300)
  set.seed(45)
  alias <- sample_posterior(post, n = 300)

  expect_equal(alias, direct)
})

test_that("sample.default still forwards to base::sample", {
  # The package masks base::sample with a generic, so the fallback path is part
  # of the contract: sample() on anything without a method must behave as before.
  set.seed(46)
  permuted <- sample(1:10)
  set.seed(46)
  expect_equal(permuted, base::sample(1:10))

  set.seed(47)
  drawn <- sample(letters, 5, replace = TRUE)
  set.seed(47)
  expect_equal(drawn, base::sample(letters, 5, replace = TRUE))

  set.seed(48)
  n_only <- sample(5)
  set.seed(48)
  expect_equal(n_only, base::sample(5))
})
