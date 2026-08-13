stan_prior <- function() prior_uniform(c(mu = -3, nu = -3), c(mu = 3, nu = 3))

stan_sim <- function(mu, nu) {
  c(a = mu + stats::rnorm(1, sd = 0.4),
    b = nu + 0.5 * mu + stats::rnorm(1, sd = 0.3))
}

stan_lingauss_fit <- function(n = 2000, seed = 1) {
  nle(stan_prior(), stan_sim, n_simulations = n,
      density_estimator = "linear_gaussian", seed = seed)
}

test_that("stan_code() emits the blocks a runnable model needs", {
  code <- stan_code(stan_lingauss_fit())

  expect_type(code, "character")
  expect_length(code, 1L)
  expect_match(code, "functions \\{")
  expect_match(code, "real nsbi_log_lik_lpdf\\(vector x, vector theta, vector w\\)")
  expect_match(code, "real nsbi_log_lik_sum_lpdf\\(matrix x, vector theta, vector w\\)")
  expect_match(code, "data \\{")
  expect_match(code, "parameters \\{")
  expect_match(code, "x ~ nsbi_log_lik_sum\\(theta, nsbi_w\\);")
})

test_that("model = FALSE gives the functions block alone", {
  fns <- stan_code(stan_lingauss_fit(), model = FALSE)

  expect_match(fns, "functions \\{")
  expect_false(grepl("parameters \\{", fns))
  expect_false(grepl("^data \\{", fns))
})

test_that("the generated names follow `name`", {
  code <- stan_code(stan_lingauss_fit(), name = "my_lik", model = FALSE)

  expect_match(code, "real my_lik_lpdf\\(")
  expect_match(code, "real my_lik_sum_lpdf\\(")
  expect_error(stan_code(stan_lingauss_fit(), name = "2bad"),
               "valid Stan identifier")
})

test_that("code generation is deterministic", {
  fit <- stan_lingauss_fit()
  expect_identical(stan_code(fit), stan_code(fit))
})

test_that("stan_data() carries the weights, the data and the prior", {
  fit <- stan_lingauss_fit()
  x_obs <- matrix(stats::rnorm(12), ncol = 2)
  data <- stan_data(fit, x_obs)

  expect_named(data, c("nsbi_nw", "nsbi_w", "N", "x", "nsbi_low", "nsbi_high"))
  expect_equal(data$nsbi_nw, length(data$nsbi_w))
  expect_equal(data$N, 6L)
  expect_equal(dim(data$x), c(6L, 2L))
  expect_equal(data$nsbi_low, c(-3, -3))
  expect_equal(data$nsbi_high, c(3, 3))
  expect_null(dimnames(data$x))
})

test_that("a normal prior is written out as a sampling statement", {
  fit <- nle(prior_normal(mean = c(mu = 0, nu = 0), sd = c(1, 2)), stan_sim,
             n_simulations = 800, density_estimator = "linear_gaussian",
             seed = 2)
  code <- stan_code(fit)

  expect_match(code, "theta ~ normal\\(nsbi_prior_mean, nsbi_prior_sd\\);")
  expect_equal(stan_data(fit)$nsbi_prior_sd, c(1, 2))
})

test_that("a custom prior cannot be written out, and the error says why", {
  prior <- prior_custom(
    sample_fn = function(n) matrix(stats::runif(2 * n, -2, 2), ncol = 2),
    log_prob_fn = function(theta) rep(-log(16), nrow(theta)),
    dim = 2L, lower = c(-2, -2), upper = c(2, 2))
  # An unnamed prior hands the simulator one vector, not named arguments.
  fit <- nle(prior, function(theta) stan_sim(theta[1], theta[2]),
             n_simulations = 500, density_estimator = "linear_gaussian",
             seed = 3)

  expect_error(stan_code(fit), "write the model block yourself")
  # The functions block is still generated: that is the point of the escape.
  expect_match(stan_code(fit, model = FALSE), "real nsbi_log_lik_lpdf\\(")
})

test_that("write_stan_model() writes a file", {
  fit <- stan_lingauss_fit()
  path <- tempfile(fileext = ".stan")
  on.exit(unlink(path), add = TRUE)

  expect_equal(write_stan_model(fit, path), path)
  expect_true(file.exists(path))
  expect_match(paste(readLines(path), collapse = "\n"), "nsbi_log_lik_sum_lpdf")
})

test_that("write_stan_model() checks `file` before it generates any code", {
  fit <- stan_lingauss_fit()

  expect_error(write_stan_model(fit, c("a", "b")),
               "`file` must be a single file path")
  expect_error(write_stan_model(fit, NULL), "`file` must be a single file path")
  expect_error(write_stan_model(fit, ""), "an empty string")
})

# cmdstan_ready() and the backend dispatch in stan_sample_nle() are what
# issue #172 is about: requireNamespace("cmdstanr") alone is TRUE whenever the
# R package is installed, even when the separate CmdStan toolchain
# (cmdstanr::install_cmdstan()) never was. These tests mock the checks rather
# than needing cmdstanr/CmdStan/rstan actually present, so they run in every
# environment, including this one, which has neither.

test_that("cmdstan_ready() is FALSE when the cmdstanr package itself is missing", {
  local_mocked_bindings(requireNamespace = function(package, ...) FALSE,
                        .package = "base")
  expect_false(cmdstan_ready())
})

test_that("cmdstan_ready() is FALSE when cmdstanr is installed but CmdStan is not", {
  # Only meaningful with cmdstanr actually present, so cmdstan_version() can be
  # mocked at all; otherwise it is covered indirectly by the test above, since
  # requireNamespace("cmdstanr") is FALSE on a machine without the package.
  testthat::skip_if_not_installed("cmdstanr")
  local_mocked_bindings(cmdstan_version = function(...) stop("path not set"),
                        .package = "cmdstanr")
  expect_false(cmdstan_ready())
})

test_that("stan_sample_nle() uses cmdstanr when cmdstan_ready() says yes", {
  fit <- stan_lingauss_fit()
  ctl <- list(n_chains = 2L, warmup = 10L, thin = 1L, seed = 1L, dots = list())
  x_obs <- matrix(stats::rnorm(4), ncol = 2)
  local_mocked_bindings(cmdstan_ready = function() TRUE)
  local_mocked_bindings(stan_run_cmdstanr = function(...) array(0, dim = c(2L, 2L, 2L)))

  out <- stan_sample_nle(fit, x_obs, ctl, n = 4)
  expect_equal(dim(out$draws), c(4L, 2L))
})

test_that("stan_sample_nle() falls back to rstan and says why, when cmdstanr isn't usable", {
  fit <- stan_lingauss_fit()
  ctl <- list(n_chains = 2L, warmup = 10L, thin = 1L, seed = 1L, dots = list())
  x_obs <- matrix(stats::rnorm(4), ncol = 2)
  local_mocked_bindings(cmdstan_ready = function() FALSE)
  local_mocked_bindings(requireNamespace = function(package, ...) identical(package, "rstan"),
                        .package = "base")
  local_mocked_bindings(stan_run_rstan = function(...) array(0, dim = c(2L, 2L, 2L)))

  expect_message(out <- stan_sample_nle(fit, x_obs, ctl, n = 4),
                 "falling back to rstan")
  expect_equal(dim(out$draws), c(4L, 2L))
})

test_that("stan_sample_nle() names install_cmdstan() when only the cmdstanr package is present", {
  fit <- stan_lingauss_fit()
  ctl <- list(n_chains = 2L, warmup = 10L, thin = 1L, seed = 1L, dots = list())
  x_obs <- matrix(stats::rnorm(4), ncol = 2)
  local_mocked_bindings(cmdstan_ready = function() FALSE)
  local_mocked_bindings(requireNamespace = function(package, ...) identical(package, "cmdstanr"),
                        .package = "base")

  expect_error(stan_sample_nle(fit, x_obs, ctl, n = 4), "install_cmdstan\\(\\)")
})

test_that("stan_sample_nle() names both packages when neither is present", {
  fit <- stan_lingauss_fit()
  ctl <- list(n_chains = 2L, warmup = 10L, thin = 1L, seed = 1L, dots = list())
  x_obs <- matrix(stats::rnorm(4), ncol = 2)
  local_mocked_bindings(cmdstan_ready = function() FALSE)
  local_mocked_bindings(requireNamespace = function(package, ...) FALSE,
                        .package = "base")

  expect_error(stan_sample_nle(fit, x_obs, ctl, n = 4),
               "needs cmdstanr or rstan installed")
})

test_that("an NSF fit is refused with the alternatives named", {
  skip_if_no_torch()
  fit <- nle(stan_prior(), stan_sim, n_simulations = 400,
             density_estimator = "nsf", hidden = c(8L, 8L),
             n_transforms = 2L, max_epochs = 5L, seed = 4)

  expect_error(stan_code(fit), "Cannot export a 'nsf' estimator")
  expect_error(stan_code(fit), "refit with \"maf\"")
})

test_that("stan_data() points a dead fit at save_npe(), as stan_code() does", {
  # readRDS() on a torch-backed fit returns a module whose external pointer is
  # nil. stan_code() has said so since it was written; stan_data() reached
  # net_param() and failed on the pointer instead.
  fit <- stan_lingauss_fit(400)
  dead <- structure(list(), class = "nsbi_dead_stan_net")
  registerS3method("$", "nsbi_dead_stan_net",
                   function(x, name) stop("external pointer is not valid"))
  fit$de$net <- dead

  expect_error(stan_data(fit), "save_npe")
  expect_error(stan_code(fit), "save_npe")
  # A live fit is untouched by the check.
  expect_type(stan_data(stan_lingauss_fit(400))$nsbi_w, "double")
})

test_that("stan_code() refuses an NPE fit", {
  fit <- npe(stan_prior(), stan_sim, n_simulations = 400,
             density_estimator = "linear_gaussian", seed = 5)
  expect_error(stan_code(fit), "needs a fit from nle\\(\\)")
})

test_that("stan_code() refuses an NRE fit and says why", {
  fit <- nre(stan_prior(), stan_sim, n_simulations = 400,
             classifier = "logistic", seed = 5)
  expect_error(stan_code(fit), "no Stan export for a ratio estimator")
  expect_error(stan_data(fit), "no Stan export for a ratio estimator")
})

test_that("the generated Stan agrees with log_lik() for linear_gaussian", {
  skip_if_no_cmdstan()
  fit <- stan_lingauss_fit()
  set.seed(6)
  theta <- matrix(stats::runif(20, -2, 2), ncol = 2)
  x_obs <- matrix(stats::rnorm(14), ncol = 2)

  from_stan <- stan_eval_log_lik(fit, theta, x_obs)

  # Closed form on both sides, so this is a pure double-precision comparison.
  expect_equal(from_stan$sum, log_lik(fit, theta, x_obs), tolerance = 1e-10)
  expect_equal(from_stan$one, log_lik(fit, theta, x_obs[1, , drop = FALSE]),
               tolerance = 1e-10)
})

test_that("the generated Stan agrees with log_lik() for MDN and MAF", {
  skip_if_no_torch()
  skip_if_no_cmdstan()
  set.seed(7)
  theta <- matrix(stats::runif(12, -2, 2), ncol = 2)
  x_obs <- matrix(stats::rnorm(10), ncol = 2)

  for (estimator in c("mdn", "maf")) {
    fit <- nle(stan_prior(), stan_sim, n_simulations = 800,
               density_estimator = estimator, hidden = c(16L, 16L),
               n_components = 3L, n_transforms = 3L, max_epochs = 12L, seed = 8)
    from_stan <- stan_eval_log_lik(fit, theta, x_obs)

    # Not exact: the network trains and evaluates in float32 while the
    # generated Stan works in double, so the two differ at single precision.
    expect_equal(from_stan$sum, log_lik(fit, theta, x_obs),
                 tolerance = 1e-4, label = paste(estimator, "sum"))
    expect_equal(from_stan$one, log_lik(fit, theta, x_obs[1, , drop = FALSE]),
                 tolerance = 1e-4, label = paste(estimator, "single"))
  }
})

test_that("stan_code() generates the MDN and MAF _sum_lpdf block without cmdstan", {
  # stan_fn_mdn()/stan_fn_maf() are only otherwise exercised by the
  # cmdstan-gated numeric round-trip above, which neither CI job runs (no
  # CmdStan installed). This checks the generated text alone, so it still
  # runs wherever torch does, including the coverage job.
  skip_if_no_torch()

  for (estimator in c("mdn", "maf")) {
    fit <- nle(stan_prior(), stan_sim, n_simulations = 400,
               density_estimator = estimator, hidden = c(8L, 8L),
               n_components = 2L, n_transforms = 2L, max_epochs = 5L, seed = 9)
    code <- stan_code(fit, model = FALSE)

    expect_match(code, "real nsbi_log_lik_sum_lpdf\\(matrix x, vector theta, vector w\\)",
                label = paste(estimator, "sum_lpdf signature"))
    expect_match(code, "for \\(n in 1:rows\\(x\\)\\) \\{", label = paste(estimator, "loop"))
    expect_match(code, "real total = 0;", label = paste(estimator, "accumulator"))
  }
})

test_that("the MAF export handles a one-dimensional observation", {
  # With dim_x = 1 the flow skips its order reversal, a branch the
  # two-dimensional cases never reach.
  skip_if_no_torch()
  skip_if_no_cmdstan()
  set.seed(9)
  prior <- prior_uniform(c(mu = -3), c(mu = 3))
  fit <- nle(prior, function(mu) c(y = stats::rnorm(1, mu, 0.5)),
             n_simulations = 800, density_estimator = "maf",
             hidden = c(16L, 16L), n_transforms = 3L, max_epochs = 12L,
             seed = 10)
  theta <- matrix(stats::runif(6, -2, 2), ncol = 1)
  x_obs <- matrix(stats::rnorm(5), ncol = 1)

  expect_equal(stan_eval_log_lik(fit, theta, x_obs)$sum,
               log_lik(fit, theta, x_obs), tolerance = 1e-4)
})

test_that("Stan NUTS and the slice sampler reach the same posterior", {
  skip_if_no_cmdstan()
  set.seed(10)
  task <- task_gaussian_linear(dim = 2L)
  fit <- nle(task$prior, task$simulator, n_simulations = 4000,
             density_estimator = "linear_gaussian", seed = 11)
  x_obs <- matrix(task$simulator(c(0.2, -0.1)), nrow = 1)

  by_stan <- sample(posterior(fit, x_obs, sampler = "stan", n_chains = 4,
                              seed = 12), 2000)
  reference <- task$reference_posterior(x_obs, 2000)

  # Absolute, not relative: these means sit near zero, where testthat's
  # relative tolerance is far tighter than the Monte-Carlo error.
  expect_lt(max(abs(colMeans(by_stan) - colMeans(reference))), 0.05)
  expect_lt(c2st(by_stan, reference, seed = 1)$accuracy, 0.6)
})

test_that("the rstan fallback compiles, samples, and hands back the same shape", {
  # The dispatch prefers cmdstanr, so the fallback is unreachable through
  # posterior() on a machine that has both. Calling it directly is the only way
  # to find out whether the branch works.
  skip_if_no_rstan()
  set.seed(10)
  task <- task_gaussian_linear(dim = 2L)
  fit <- nle(task$prior, task$simulator, n_simulations = 4000,
             density_estimator = "linear_gaussian", seed = 11)
  x_obs <- matrix(task$simulator(c(0.2, -0.1)), nrow = 1)

  draws <- stan_run_rstan(stan_code(fit), stan_data(fit, x_obs),
                          ctl = list(n_chains = 2L, seed = 12),
                          iter_warmup = 400L, iter_sampling = 400L,
                          refresh = 0L)

  # iterations x chains x dim, which is what mcmc_diagnostics() reads and what
  # the cmdstanr branch produces.
  expect_equal(dim(draws), c(400L, 2L, 2L))
  expect_true(all(is.finite(draws)))
  expect_true(all(mcmc_diagnostics(draws)$rhat < 1.05))

  flat <- matrix(aperm(draws, c(2, 1, 3)), ncol = 2)
  reference <- task$reference_posterior(x_obs, 2000)
  expect_lt(max(abs(colMeans(flat) - colMeans(reference))), 0.06)
})
