test_that("an NLE posterior matches the analytic one on the conjugate task", {
  # x | theta is Gaussian and linear in theta, so the linear_gaussian
  # estimator is exact and any disagreement is the sampler's fault.
  set.seed(1)
  task <- task_gaussian_linear(dim = 3L)
  fit <- nle(task$prior, task$simulator, n_simulations = 4000,
             density_estimator = "linear_gaussian", seed = 1)
  x_obs <- matrix(task$simulator(c(0.1, -0.2, 0.3)), nrow = 1)

  post <- posterior(fit, x_obs, n_chains = 10, warmup = 100, thin = 5, seed = 2)
  draws <- sample(post, 3000)
  reference <- task$reference_posterior(x_obs, 3000)

  expect_s3_class(post, "nsbi_nle_posterior")
  expect_s3_class(post, "nsbi_posterior")
  expect_s3_class(draws, "nsbi_samples")
  # Absolute, not relative: these means sit near zero, where testthat's
  # relative tolerance is far tighter than the Monte-Carlo error.
  expect_lt(max(abs(colMeans(draws) - colMeans(reference))), 0.05)
  expect_equal(apply(draws, 2, stats::sd), apply(reference, 2, stats::sd),
               tolerance = 0.1)
  expect_lt(c2st(draws, reference, seed = 1)$accuracy, 0.6)
})

test_that("more independent observations tighten the posterior correctly", {
  # The reason to prefer NLE here: one fit, any number of trials. The
  # conjugate answer for n observations has variance 1/(1/v0 + n/s2).
  set.seed(2)
  prior_var <- 0.1
  noise_var <- 0.1
  prior <- prior_normal(mean = 0, sd = sqrt(prior_var))
  simulator <- function(theta) stats::rnorm(1, theta, sd = sqrt(noise_var))
  fit <- nle(prior, simulator, n_simulations = 5000,
             density_estimator = "linear_gaussian", seed = 3)

  for (n_obs in c(1L, 10L, 100L)) {
    x_obs <- matrix(stats::rnorm(n_obs, mean = 0.3, sd = sqrt(noise_var)), ncol = 1)
    post <- posterior(fit, x_obs, n_chains = 8, warmup = 100, thin = 4, seed = 4)
    draws <- sample(post, 2000)

    post_var <- 1 / (1 / prior_var + n_obs / noise_var)
    post_mean <- post_var * sum(x_obs) / noise_var

    expect_lt(abs(mean(draws) - post_mean), 0.05,
              label = sprintf("mean with n_obs = %d", n_obs))
    expect_equal(stats::sd(as.numeric(draws)), sqrt(post_var), tolerance = 0.1,
                 label = sprintf("sd with n_obs = %d", n_obs))
  }
})

test_that("draws are cached, and a seeded posterior is reproducible", {
  set.seed(3)
  prior <- prior_uniform(c(mu = -3), c(mu = 3))
  fit <- nle(prior, function(mu) c(y = stats::rnorm(1, mu, 0.5)),
             n_simulations = 1500, density_estimator = "linear_gaussian",
             seed = 5)
  post <- posterior(fit, matrix(stats::rnorm(20, 1, 0.5), ncol = 1),
                    n_chains = 4, warmup = 50, thin = 2, seed = 6)

  first <- sample(post, 400)
  again <- sample(post, 400)
  fewer <- sample(post, 200)

  expect_equal(again, first, ignore_attr = TRUE)
  expect_equal(fewer, first[1:200, , drop = FALSE], ignore_attr = TRUE)
  # The seed is re-applied per run, so refresh reproduces rather than redraws.
  expect_equal(sample(post, 400, refresh = TRUE), first, ignore_attr = TRUE)
})

test_that("refresh re-runs the chain when the posterior is unseeded", {
  set.seed(30)
  prior <- prior_uniform(c(mu = -3), c(mu = 3))
  fit <- nle(prior, function(mu) c(y = stats::rnorm(1, mu, 0.5)),
             n_simulations = 1500, density_estimator = "linear_gaussian",
             seed = 5)
  post <- posterior(fit, matrix(stats::rnorm(20, 1, 0.5), ncol = 1),
                    n_chains = 4, warmup = 50, thin = 2)

  first <- sample(post, 400)

  expect_equal(sample(post, 400), first, ignore_attr = TRUE)
  expect_false(isTRUE(all.equal(unname(unclass(sample(post, 400, refresh = TRUE))),
                                unname(unclass(first)))))
})

test_that("a new observation is not served from the cache", {
  set.seed(4)
  prior <- prior_uniform(c(mu = -3), c(mu = 3))
  fit <- nle(prior, function(mu) c(y = stats::rnorm(1, mu, 0.5)),
             n_simulations = 2000, density_estimator = "linear_gaussian",
             seed = 7)
  post <- posterior(fit, n_chains = 6, warmup = 50, thin = 2, seed = 8)

  low <- sample(post, 600, obs = matrix(stats::rnorm(50, -1, 0.5), ncol = 1))
  high <- sample(post, 600, obs = matrix(stats::rnorm(50, 1.5, 0.5), ncol = 1))

  expect_lt(mean(low), mean(high) - 1)
})

test_that("posterior draws carry names and convergence diagnostics", {
  set.seed(5)
  prior <- prior_uniform(c(mu = -3, nu = -3), c(mu = 3, nu = 3))
  fit <- nle(prior,
             function(mu, nu) c(a = mu + stats::rnorm(1, sd = 0.4),
                                b = nu + stats::rnorm(1, sd = 0.4)),
             n_simulations = 2000, density_estimator = "linear_gaussian",
             seed = 9)
  post <- posterior(fit, matrix(c(0.5, -0.5), nrow = 1),
                    n_chains = 6, warmup = 60, thin = 3, seed = 10)
  draws <- sample(post, 900)
  diagnostics <- attr(draws, "diagnostics")

  expect_equal(colnames(draws), c("mu", "nu"))
  expect_equal(rownames(diagnostics), c("mu", "nu"))
  expect_true(all(diagnostics$rhat < 1.05))
  expect_output(print(post), "sampler         : slice")
  expect_output(print(post), "unnormalized")
})

test_that("print() says when the last run has no usable diagnostics", {
  set.seed(6)
  prior <- prior_uniform(c(mu = -3), c(mu = 3))
  fit <- nle(prior, function(mu) c(y = stats::rnorm(1, mu, 0.5)),
             n_simulations = 800, density_estimator = "linear_gaussian",
             seed = 11)
  post <- posterior(fit, matrix(0.5, nrow = 1), n_chains = 2, warmup = 20,
                    seed = 12)
  sample(post, 100)

  # What mcmc_diagnostics() returns for a run too short or too degenerate to
  # score. max(na.rm = TRUE) over it is -Inf, which is what used to print.
  post$cache$diagnostics <- data.frame(rhat = NA_real_, ess_bulk = NA_real_)
  out <- utils::capture.output(print(post))
  expect_true(any(grepl("diagnostics unavailable", out)))
  expect_false(any(grepl("Inf", out, fixed = TRUE)))
})

test_that("summary() works through the sample generic", {
  set.seed(20)
  prior <- prior_uniform(c(mu = -3), c(mu = 3))
  fit <- nle(prior, function(mu) c(y = stats::rnorm(1, mu, 0.5)),
             n_simulations = 1500, density_estimator = "linear_gaussian",
             seed = 21)
  x_obs <- matrix(stats::rnorm(30, 1, 0.5), ncol = 1)
  post <- posterior(fit, x_obs, n_chains = 4, warmup = 50, seed = 22)

  s <- summary(post, n = 500)

  expect_s3_class(s, "data.frame")
  expect_equal(s$parameter, "mu")
  expect_lt(abs(s$mean - mean(x_obs)), 0.15)
})

test_that("map_estimate() finds the mode of an NLE posterior", {
  # map_estimate() used to call sample.nsbi_posterior() rather than the
  # generic, which on an nle() fit asks the estimator for draws with the roles
  # swapped: it returns points in x space. With dim_x = 1 and dim_theta = 2
  # that is a non-conformable error; where the two dimensions happen to match
  # it is silently the wrong starting distribution.
  set.seed(31)
  prior <- prior_uniform(c(mu = -3, nu = -3), c(mu = 3, nu = 3))
  simulator <- function(mu, nu) c(y = mu + 0.5 * nu + stats::rnorm(1, sd = 0.3))
  fit <- nle(prior, simulator, n_simulations = 3000,
             density_estimator = "linear_gaussian", seed = 32)
  x_obs <- matrix(stats::rnorm(50, mean = 1, sd = 0.3), ncol = 1)
  post <- posterior(fit, x_obs, n_chains = 6, warmup = 60, thin = 2, seed = 33)

  map <- map_estimate(post, n_init = 400)

  expect_length(map, 2L)
  expect_equal(names(map), c("mu", "nu"))
  # Only mu + nu/2 is identified, so the ridge is what can be checked.
  expect_equal(unname(map[1] + 0.5 * map[2]), mean(x_obs), tolerance = 0.1)
  # And the mode must beat the posterior mean at its own game.
  expect_gte(log_prob(post, matrix(map, nrow = 1), normalize = FALSE),
             log_prob(post, matrix(colMeans(sample(post, 400)), nrow = 1),
                      normalize = FALSE) - 1e-6)
})

test_that("log_prob is the unnormalized posterior and says so", {
  set.seed(6)
  prior <- prior_uniform(c(mu = -2), c(mu = 2))
  fit <- nle(prior, function(mu) c(y = stats::rnorm(1, mu, 0.5)),
             n_simulations = 1500, density_estimator = "linear_gaussian",
             seed = 11)
  x_obs <- matrix(stats::rnorm(10, 0.5, 0.5), ncol = 1)
  post <- posterior(fit, x_obs, n_chains = 4)

  theta <- matrix(c(0.5, 1.5, 5), ncol = 1)
  lp <- log_prob(post, theta, normalize = FALSE)

  expect_length(lp, 3L)
  expect_equal(lp[1:2],
               log_lik(fit, theta[1:2, , drop = FALSE], x_obs) +
                 prior$log_prob(theta[1:2, , drop = FALSE]))
  expect_equal(lp[3], -Inf)               # outside the prior support
  expect_warning(log_prob(post, theta, normalize = TRUE), "unnormalized")
})

test_that("sampling refuses a prior with no density", {
  set.seed(7)
  prior <- prior_custom(sample_fn = function(n) matrix(stats::runif(n, -2, 2), ncol = 1),
                        dim = 1L, lower = -2, upper = 2)
  fit <- nle(prior, function(theta) stats::rnorm(1, theta, 0.5),
             n_simulations = 800, density_estimator = "linear_gaussian",
             seed = 12)
  post <- posterior(fit, matrix(0, nrow = 1), n_chains = 4)

  expect_error(sample(post, 100), "log_prob_fn")
})

test_that("a single chain is refused, since it cannot be diagnosed", {
  prior <- prior_uniform(c(mu = -2), c(mu = 2))
  fit <- nle(prior, function(mu) c(y = stats::rnorm(1, mu, 0.5)),
             n_simulations = 500, density_estimator = "linear_gaussian")
  expect_error(posterior(fit, matrix(0), n_chains = 1), "at least 2")
})

test_that("thin and warmup are validated, and thin = 1 keeps every draw", {
  prior <- prior_uniform(c(mu = -2), c(mu = 2))
  fit <- nle(prior, function(mu) c(y = stats::rnorm(1, mu, 0.5)),
             n_simulations = 500, density_estimator = "linear_gaussian",
             seed = 40)
  x_obs <- matrix(stats::rnorm(10, 0.5, 0.5), ncol = 1)

  # thin = 0 used to run warmup only and hand back the zero-initialized array
  # of kept draws, with no error and a plausible-looking rhat.
  expect_error(posterior(fit, x_obs, thin = 0), "`thin` must be")
  expect_error(posterior(fit, x_obs, thin = -1), "`thin` must be")
  expect_error(posterior(fit, x_obs, thin = 1.5), "`thin` must be")
  expect_error(posterior(fit, x_obs, thin = NA), "`thin` must be")
  expect_error(posterior(fit, x_obs, warmup = -5), "`warmup` must be")

  post <- posterior(fit, x_obs, n_chains = 4, warmup = 40, thin = 1, seed = 41)
  draws <- sample(post, 200)
  expect_equal(nrow(draws), 200L)
  expect_false(all(draws == 0))
  expect_gt(stats::sd(as.numeric(draws)), 0)
})

test_that("posterior() rejects an object that is neither fit", {
  expect_error(posterior(list(a = 1)), "needs a fit from npe\\(\\) or nle\\(\\)")
})

test_that("sbc() and posterior_predictive() accept an NLE fit", {
  set.seed(8)
  prior <- prior_uniform(c(mu = -2), c(mu = 2))
  simulator <- function(mu) c(y = stats::rnorm(1, mu, 0.5))
  fit <- nle(prior, simulator, n_simulations = 3000,
             density_estimator = "linear_gaussian", seed = 13)

  res <- sbc(fit, simulator, n_sbc = 30, n_posterior_samples = 100, seed = 14,
             n_chains = 4, warmup = 40, thin = 2)
  expect_s3_class(res, "nsbi_sbc")
  expect_equal(nrow(res$ranks), 30L)
  # A correct fit gives uniform ranks; 30 trials only catches gross failure.
  expect_gt(res$uniformity_pvalue[[1]], 0.001)

  post <- posterior(fit, matrix(0.5, nrow = 1), n_chains = 4, warmup = 40,
                    thin = 2, seed = 15)
  pred <- posterior_predictive(post, simulator, n = 100)
  expect_equal(dim(pred), c(100L, 1L))
})

test_that("a neural NLE posterior lands near the analytic one", {
  skip_if_no_torch()
  set.seed(9)
  task <- task_gaussian_linear(dim = 2L)
  fit <- nle(task$prior, task$simulator, n_simulations = 5000,
             density_estimator = "maf", max_epochs = 300L, seed = 16)
  x_obs <- matrix(task$simulator(c(0.2, -0.1)), nrow = 1)

  post <- posterior(fit, x_obs, n_chains = 10, warmup = 150, thin = 5, seed = 17)
  draws <- sample(post, 2000)
  reference <- task$reference_posterior(x_obs, 2000)

  # Looser than the closed-form case: a trained flow is an approximation.
  expect_lt(max(abs(colMeans(draws) - colMeans(reference))), 0.1)
  expect_equal(apply(draws, 2, stats::sd), apply(reference, 2, stats::sd),
               tolerance = 0.25)
})
