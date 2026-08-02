# The simulator contract: how one parameter set reaches the simulator, and
# what it is allowed to give back. See ?nsbi_simulator.

test_that("parameters arrive by name when they match the formals", {
  prior <- prior_uniform(low = c(mu = 0, sigma = 1),
                         high = c(mu = 1, sigma = 2))
  expect_equal(sim_dispatch(function(mu, sigma) NULL, c("mu", "sigma")),
               "named")
  # a spare formal is fine; a missing one is not
  expect_equal(sim_dispatch(function(mu, sigma, n) NULL, c("mu", "sigma")),
               "named")
  expect_warning(
    expect_equal(sim_dispatch(function(mu) NULL, c("mu", "sigma")), "vector"),
    "sigma has no formal")
  expect_equal(sim_dispatch(function(theta) NULL, c("mu", "sigma")), "vector")
  # `...` cannot stand in for a parameter: the match has to be explicit
  expect_equal(sim_dispatch(function(...) NULL, c("mu", "sigma")), "vector")
  # an unnamed prior always takes the vector form
  expect_equal(sim_dispatch(function(mu, sigma) NULL, NULL), "vector")

  sims <- simulate_for_sbi(function(mu, sigma) c(m = mu, s = sigma),
                           prior, 20, seed = 1)
  expect_equal(colnames(sims$x), c("m", "s"))
  expect_equal(unname(sims$x), unname(sims$theta))
})

test_that("a partial match between parameter names and formals warns", {
  # Two of three formals match, so the vector form sends the whole parameter
  # vector to `mu` and `ls` keeps its default. Silently fitting on that is the
  # bug this warning exists for.
  expect_warning(
    expect_equal(
      sim_dispatch(function(mu, sigma, ls = 0) NULL, c("mu", "sigma", "tau")),
      "vector"),
    "tau has no formal")
  expect_warning(
    sim_dispatch(function(mu, ls = 0) NULL, c("a", "mu")),
    "Passing the whole parameter vector to `mu`", fixed = TRUE)
  expect_warning(
    sim_dispatch(function(mu, ls = 0) NULL, c("mu", "a", "b")),
    "a, b have no formal")

  prior <- prior_uniform(low = c(a = -3, b = -1), high = c(a = 3, b = 1))
  expect_warning(
    simulate_for_sbi(function(a, ls = 0) rnorm(1, a, exp(ls)), prior, 5,
                     seed = 1),
    "b has no formal")

  # a plain vector-signature simulator stays silent
  expect_silent(sim_dispatch(function(theta) NULL, c("mu", "sigma")))
  expect_silent(sim_dispatch(function(mu, sigma) NULL, c("mu", "sigma")))
  expect_silent(sim_dispatch(function(...) NULL, c("mu", "sigma")))
})

test_that("a non-syntactic parameter name always takes the vector form", {
  prior <- prior_normal(mean = c(`beta[1]` = 0, rho = 0), sd = 1)
  expect_equal(sim_dispatch(function(theta) NULL, prior$param_names), "vector")
  sims <- simulate_for_sbi(function(theta) unname(theta), prior, 20, seed = 1)
  expect_equal(unname(sims$x), unname(sims$theta))
})

test_that("sim_args reach the simulator under both signatures", {
  prior <- prior_uniform(low = c(alpha = -1, beta = -1),
                         high = c(alpha = 1, beta = 1))
  grid <- c(0, 0.5, 1)

  by_name <- function(alpha, beta, x_grid) alpha + beta * x_grid
  x <- simulate_for_sbi(by_name, prior, 10, sim_args = list(x_grid = grid),
                        seed = 2)
  expect_equal(dim(x$x), c(10L, 3L))
  expect_equal(x$x[, 1], unname(x$theta[, "alpha"]))

  by_vector <- function(theta, x_grid) theta["alpha"] + theta["beta"] * x_grid
  y <- simulate_for_sbi(by_vector, prior, 10, sim_args = list(x_grid = grid),
                        seed = 2)
  expect_equal(x$x, y$x)
})

test_that("sim_args is checked before anything runs", {
  prior <- prior_uniform(low = c(a = 0), high = c(a = 1))
  expect_error(
    simulate_for_sbi(function(a, b) a + b, prior, 5, sim_args = list(1)),
    "must be named"
  )
  expect_error(
    simulate_for_sbi(function(a) a, prior, 5, sim_args = list(a = 1)),
    "clash with parameter names: a"
  )
  expect_error(simulate_for_sbi("not a function", prior, 5),
               "must be a function")
})

test_that("a single named parameter keeps its name -- the documented sharp edge", {
  prior <- prior_uniform(low = c(beta = 1, gamma = 1),
                         high = c(beta = 2, gamma = 2))
  # theta["beta"] / theta["gamma"] is a scalar still called "beta"
  sims <- simulate_for_sbi(function(theta) theta["beta"] / theta["gamma"],
                           prior, 10, seed = 1)
  expect_equal(colnames(sims$x), "beta")
  # unname() is the fix
  sims <- simulate_for_sbi(
    function(theta) unname(theta["beta"] / theta["gamma"]), prior, 10, seed = 1)
  expect_null(colnames(sims$x))
})

test_that("a failure inside the simulator names the simulation", {
  prior <- prior_uniform(low = c(mu = -3, sigma = 0.1),
                         high = c(mu = 3, sigma = 3))
  calls <- 0L
  sim <- function(mu, sigma) {
    calls <<- calls + 1L
    if (calls == 3L) stop("the solver did not converge")
    c(y = mu)
  }
  err <- expect_error(simulate_for_sbi(sim, prior, 5, seed = 1))
  expect_match(conditionMessage(err),
               "Simulation 3 failed: the solver did not converge", fixed = TRUE)
  # the parameters that produced it, by name, rounded
  expect_match(conditionMessage(err), "parameters: mu = ", fixed = TRUE)
  expect_match(conditionMessage(err), "sigma = ", fixed = TRUE)
  expect_match(conditionMessage(err), "?nsbi_simulator", fixed = TRUE)
})

test_that("a dispatch mismatch is reported as one, not as a missing default", {
  # unnamed prior, so the whole vector goes to `mu` and `ls` has no default.
  # The bare R error names neither the simulation nor the contract.
  prior <- prior_uniform(c(-3, -1), c(3, 1))
  err <- expect_error(
    simulate_for_sbi(function(mu, ls) rnorm(1, mu, exp(ls)), prior, 5, seed = 1))
  expect_match(conditionMessage(err), "Simulation 1 failed", fixed = TRUE)
  expect_match(conditionMessage(err), "?nsbi_simulator", fixed = TRUE)
  # no parameter names to print, so the values stand on their own
  expect_match(conditionMessage(err), "parameters: -?[0-9]")
})

test_that("the simulator's own condition survives the re-raise", {
  cnd <- errorCondition("diverged", class = "my_solver_error")
  call_one <- function(theta_i) stop(cnd)
  err <- expect_error(call_sim_once(call_one, c(mu = 1), 2L),
                      class = "my_solver_error")
  expect_s3_class(err, "nsbi_sim_error")
  expect_identical(conditionMessage(err$parent), "diverged")
  # a call that works is passed through untouched
  expect_equal(call_sim_once(function(theta_i) theta_i * 2, c(mu = 1.5), 1L),
               c(mu = 3))
})

test_that("a wide parameter set is truncated in the message", {
  wide <- stats::setNames(seq_len(20) / 7, paste0("p", 1:20))
  err <- expect_error(call_sim_once(function(theta_i) stop("nope"), wide, 4L))
  expect_match(conditionMessage(err), "p1 = 0.1429", fixed = TRUE)
  expect_match(conditionMessage(err), "(20 parameters in all)", fixed = TRUE)
  expect_false(grepl("p20", conditionMessage(err), fixed = TRUE))
})

test_that("npe() runs end to end on a per-parameter-set simulator", {
  set.seed(11)
  prior <- prior_normal(mean = c(alpha = 0, beta = 0), sd = 1)
  simulator <- function(alpha, beta, x_grid) {
    rnorm(length(x_grid), mean = alpha + beta * x_grid, sd = 0.2)
  }
  grid <- seq(-1, 1, length.out = 4)
  fit <- npe(prior, simulator, n_simulations = 2000,
             sim_args = list(x_grid = grid),
             density_estimator = "linear_gaussian")
  expect_equal(fit$dim_x, 4L)
  expect_equal(fit$n_dropped, 0L)

  truth <- c(alpha = 0.5, beta = -0.8)
  x_obs <- truth["alpha"] + truth["beta"] * grid
  draws <- sample(posterior(fit, x_obs = x_obs), 4000)
  expect_equal(unname(colMeans(draws)), unname(truth), tolerance = 0.15)
})
