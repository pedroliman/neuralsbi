# A fit stores the column order of the simulator output as x_names. A named
# x_obs is matched to it by name, so the order the user typed does not matter
# (#389). Unnamed input is still read by position.

named_fit <- function() {
  set.seed(1)
  pr <- prior_normal(c(a = 0, b = 0), 1)
  sim <- function(a, b) {
    c(y1 = a + stats::rnorm(1, sd = 0.2), y2 = 3 * b + stats::rnorm(1, sd = 0.2))
  }
  npe(pr, sim, 600, density_estimator = "linear_gaussian")
}

test_that("a named x_obs is reordered to fit$x_names for NPE", {
  fit <- named_fit()
  expect_equal(fit$x_names, c("y1", "y2"))
  ref <- posterior(fit, x_obs = c(y1 = 1, y2 = 0))$x_obs

  vec <- posterior(fit, x_obs = c(y2 = 0, y1 = 1))$x_obs
  df <- posterior(fit, x_obs = data.frame(y2 = 0, y1 = 1))$x_obs
  mat <- posterior(fit, x_obs = matrix(c(0, 1), 1,
                                       dimnames = list(NULL, c("y2", "y1"))))$x_obs
  expect_equal(unname(vec), unname(ref))
  expect_equal(unname(df), unname(ref))
  expect_equal(unname(mat), unname(ref))
  expect_equal(colnames(vec), c("y1", "y2"))

  # the same reordering applies when the observation arrives through sample()
  post <- posterior(fit)
  set.seed(2)
  d1 <- sample(post, 500, obs = c(y1 = 1, y2 = 0))
  set.seed(2)
  d2 <- sample(post, 500, obs = c(y2 = 0, y1 = 1))
  expect_equal(unname(as.matrix(d1)), unname(as.matrix(d2)))
})

test_that("names that do not match x_names are an error", {
  fit <- named_fit()
  expect_error(posterior(fit, x_obs = c(y1 = 1, z = 0)),
               "do not match the observation names")
  expect_error(sample(posterior(fit), 10, obs = c(q = 1, y2 = 0)),
               "do not match")
})

test_that("an unnamed x_obs is still read by position", {
  fit <- named_fit()
  expect_equal(unname(posterior(fit, x_obs = c(2, 5))$x_obs),
               matrix(c(2, 5), 1))
})

test_that("NLE and NRE posteriors reorder a named x_obs too", {
  set.seed(1)
  pr <- prior_normal(c(a = 0, b = 0), 1)
  sim <- function(a, b) {
    c(y1 = a + stats::rnorm(1, sd = 0.2), y2 = 3 * b + stats::rnorm(1, sd = 0.2))
  }
  fit <- nle(pr, sim, 600, density_estimator = "linear_gaussian")
  x_obs <- rbind(c(y2 = 0, y1 = 1), c(y2 = 0.3, y1 = 0.8))
  post <- posterior(fit, x_obs = x_obs, n_chains = 2, warmup = 10)
  expect_equal(colnames(post$x_obs), c("y1", "y2"))
  expect_equal(unname(post$x_obs), unname(x_obs[, 2:1]))
  expect_error(posterior(fit, x_obs = c(y1 = 1, z = 0)), "do not match")
})

test_that("NRE posteriors reorder a named x_obs", {
  skip_if_no_torch()
  set.seed(1)
  pr <- prior_normal(c(a = 0, b = 0), 1)
  sim <- function(a, b) {
    c(y1 = a + stats::rnorm(1, sd = 0.2), y2 = 3 * b + stats::rnorm(1, sd = 0.2))
  }
  fit <- nre(pr, sim, 300, max_epochs = 2)
  post <- posterior(fit, x_obs = c(y2 = 0, y1 = 1), n_chains = 2, warmup = 10)
  expect_equal(unname(post$x_obs), matrix(c(1, 0), 1))
})

test_that("an unknown sampler argument to an NLE posterior is an error (#390)", {
  set.seed(1)
  fit <- nle(prior_uniform(-3, 3), function(mu) stats::rnorm(1, mu, .5), 600,
             density_estimator = "linear_gaussian")
  x <- matrix(1, 1, 1)
  expect_error(posterior(fit, x, n_chains = 4, wamrup = 10), "wamrup")
  expect_error(posterior(fit, x, iter_warmup = 10), "iter_warmup")
  expect_s3_class(posterior(fit, x, n_chains = 4, width = 1, n_pool = 50),
                  "nsbi_nle_posterior")
})
