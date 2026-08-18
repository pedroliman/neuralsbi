# TSNPE with the exact linear_gaussian estimator: the truncation machinery is
# exercised end to end, and the final posterior must still match the analytic
# conjugate posterior at the targeted observation.

test_that("npe_sequential recovers the analytic posterior at x_obs", {
  set.seed(21)
  d <- 2; sigma <- 0.5
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + rnorm(length(theta), sd = sigma)
  x_obs <- c(1.0, -0.5)
  fit <- npe_sequential(prior, simulator, x_obs = x_obs,
                        n_rounds = 3, n_simulations = 800,
                        density_estimator = "linear_gaussian", seed = 1)
  expect_s3_class(fit, "nsbi_snpe")
  expect_s3_class(fit, "nsbi_npe")
  expect_length(fit$rounds, 3L)
  expect_equal(fit$n_simulations, 2400L)

  post <- posterior(fit, x_obs = x_obs)
  draws <- sample(post, 5000)
  truth <- analytic_gauss_posterior(x_obs, sigma, d)
  expect_equal(colMeans(draws), truth$mu, tolerance = 0.06)
  expect_equal(apply(draws, 2, sd), sqrt(diag(truth$Sigma)), tolerance = 0.06)

  z <- matrix(rnorm(5000 * d), ncol = d)
  analytic_draws <- sweep(z %*% chol(truth$Sigma), 2, truth$mu, `+`)
  expect_lt(c2st(draws, analytic_draws, seed = 2)$accuracy, 0.6)
})

test_that("later rounds actually truncate the proposal", {
  set.seed(22)
  prior <- prior_normal(mean = 0, sd = 2)
  simulator <- function(theta) theta + rnorm(length(theta), sd = 0.2)
  fit <- npe_sequential(prior, simulator, x_obs = 1,
                        n_rounds = 2, n_simulations = 500,
                        density_estimator = "linear_gaussian", seed = 3)
  # round 1 draws from the prior itself
  expect_equal(fit$rounds[[1]]$acceptance, 1)
  expect_identical(fit$rounds[[1]]$threshold, -Inf)
  # posterior sd ~0.2 vs prior sd 2: most prior candidates must be rejected
  expect_lt(fit$rounds[[2]]$acceptance, 0.8)
  expect_gt(fit$rounds[[2]]$acceptance, 0)
  expect_true(is.finite(fit$rounds[[2]]$threshold))
  expect_output(print(fit), "nsbi_snpe")
})

test_that("per-round budgets can differ and are recorded", {
  set.seed(23)
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + rnorm(length(theta), sd = 0.5)
  fit <- npe_sequential(prior, simulator, x_obs = 0.3,
                        n_rounds = 2, n_simulations = c(600, 300),
                        density_estimator = "linear_gaussian", seed = 4)
  expect_equal(fit$rounds[[1]]$n_new, 600L)
  expect_equal(fit$rounds[[2]]$n_new, 300L)
  expect_equal(fit$n_simulations, 900L)
})

test_that("npe_sequential requires a simulator function", {
  prior <- prior_normal(0, 1)
  expect_error(npe_sequential(prior, simulator = NULL, x_obs = 0),
               "simulator")
})

test_that("npe_sequential requires x_obs to be supplied at all, not just non-NULL", {
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + rnorm(1, sd = 0.5)
  expect_error(npe_sequential(prior, simulator, n_rounds = 1, n_simulations = 100,
                              density_estimator = "linear_gaussian"),
               "`x_obs` is required")
})

test_that("npe_sequential accepts x_obs as a one-row data frame", {
  set.seed(41)
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + rnorm(1, sd = 0.5)
  fit <- npe_sequential(prior, simulator, x_obs = data.frame(y = 0.3),
                        n_rounds = 1, n_simulations = 200,
                        density_estimator = "linear_gaussian", seed = 1)
  expect_equal(fit$x_obs, 0.3)
})

test_that("npe_sequential rejects a prior that is not an nsbi_prior", {
  expect_error(npe_sequential(list(), function(theta) theta, x_obs = 0),
               "inherits")
})

test_that("npe_sequential warns and continues with fewer draws when the proposal batch cap bites", {
  set.seed(31)
  # a wide prior and a tight likelihood: round 2's truncated region is narrow,
  # so with only one proposal batch allowed, most rounds fall short of budget
  prior <- prior_uniform(low = -50, high = 50)
  simulator <- function(theta) theta + rnorm(1, sd = 0.05)
  expect_warning(
    fit <- npe_sequential(prior, simulator, x_obs = 0.2, n_rounds = 2,
                          n_simulations = 500, epsilon = 0.5,
                          max_proposal_batches = 1,
                          density_estimator = "linear_gaussian", seed = 8),
    "proposal draws inside the truncated region"
  )
  expect_lt(fit$rounds[[2]]$n_new, 500L)
})

test_that("npe_sequential errors clearly, rather than crashing rbind(), when a round accepts zero proposals and dim_x > 1", {
  set.seed(1)
  prior <- prior_uniform(low = c(-5, -5), high = c(5, 5))
  simulator <- function(theta) {
    x <- theta + rnorm(length(theta), sd = 0.01)
    matrix(x, nrow = 1)
  }
  x_obs <- matrix(c(0, 0), nrow = 1)
  # epsilon this tight plus a single proposal batch means round 2 almost
  # certainly accepts nothing; before the fix this crashed inside rbind()
  # with an opaque "number of columns of matrices must match" error because
  # run_simulator()'s zero-row fast path fabricated a 0 x 1 matrix instead
  # of a 0 x 2 one (issue #178).
  expect_error(
    suppressWarnings(npe_sequential(prior, simulator, x_obs,
                                    n_rounds = 2, n_simulations = 30,
                                    density_estimator = "linear_gaussian",
                                    epsilon = 1e-4, max_proposal_batches = 1)),
    "round 2 accepted 0/30 proposal draws"
  )
})

test_that("print.nsbi_snpe() labels the targeted x_obs by name when the fit has names", {
  set.seed(9)
  prior <- prior_normal(mean = c(beta = 0, rho = 0), sd = 1)
  simulator <- function(beta, rho) c(cases = beta, deaths = rho) + rnorm(2, sd = 0.3)
  fit <- npe_sequential(prior, simulator, x_obs = c(cases = 0.2, deaths = -0.1),
                        n_rounds = 2, n_simulations = 400,
                        density_estimator = "linear_gaussian", seed = 8)
  expect_equal(fit$x_names, c("cases", "deaths"))
  expect_output(print(fit), "cases=.*deaths=")
})

test_that("npe_sequential rejects an x_obs the simulator cannot have produced", {
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + rnorm(1, sd = 0.5)
  # dim_x is 1: a length-3 x_obs used to be folded into three rows, of which
  # only the first was ever targeted
  expect_error(
    npe_sequential(prior, simulator, x_obs = c(1, 2, 3), n_rounds = 1,
                   n_simulations = 100, density_estimator = "linear_gaussian",
                   seed = 6),
    "3 value\\(s\\) but the simulator returns 1")
  expect_error(
    npe_sequential(prior, simulator, x_obs = matrix(c(1, 2), ncol = 1),
                   n_rounds = 1, n_simulations = 100,
                   density_estimator = "linear_gaussian", seed = 6),
    "2 rows")
  expect_error(
    npe_sequential(prior, simulator, x_obs = NULL, n_rounds = 1,
                   n_simulations = 100, density_estimator = "linear_gaussian"),
    "`x_obs` is required")
  expect_error(
    npe_sequential(prior, simulator, x_obs = NA_real_, n_rounds = 1,
                   n_simulations = 100, density_estimator = "linear_gaussian"),
    "missing values")
})

test_that("npe_sequential requires at least one round", {
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + rnorm(1, sd = 0.5)
  expect_error(
    npe_sequential(prior, simulator, x_obs = 0, n_rounds = 0,
                   n_simulations = 100, density_estimator = "linear_gaussian"),
    "n_rounds")
})

test_that("truncation works with a bounded prior", {
  set.seed(24)
  prior <- prior_uniform(c(-2, -2), c(2, 2))
  simulator <- function(theta) theta + rnorm(length(theta), sd = 0.3)
  x_obs <- c(0.5, 0.5)
  fit <- npe_sequential(prior, simulator, x_obs = x_obs,
                        n_rounds = 2, n_simulations = 600,
                        density_estimator = "linear_gaussian", seed = 5)
  post <- posterior(fit, x_obs = x_obs)
  draws <- sample(post, 2000)
  # posterior concentrated near the truth, well inside the box
  expect_equal(colMeans(draws), x_obs, tolerance = 0.1)
  expect_true(all(within_support(prior, draws)))
})

test_that("npe_sequential's proposal loop requests the shortfall each batch, not the full round budget", {
  set.seed(60)
  requested <- integer(0)
  real_sample_prior <- sample_prior
  local_mocked_bindings(
    sample_prior = function(prior, n) {
      requested[[length(requested) + 1L]] <<- n
      real_sample_prior(prior, n)
    }
  )
  # a prior much wider than the likelihood, so round 2's truncated proposal
  # region is narrow enough to need several batches, not exhausting the
  # default max_proposal_batches
  prior <- prior_normal(mean = 0, sd = 3)
  simulator <- function(theta) theta + rnorm(1, sd = 0.3)
  fit <- npe_sequential(prior, simulator, x_obs = 1.5, n_rounds = 2,
                        n_simulations = 300, epsilon = 0.3,
                        density_estimator = "linear_gaussian", seed = 61)
  expect_equal(fit$rounds[[2]]$n_new, 300L)

  # requested[1] is round 1's single call (no truncation, draws straight from
  # the prior); requested[2] is round 2's first batch, which still has to ask
  # for the full budget since nothing has been collected yet
  expect_gt(length(requested), 2L)   # round 2 needed more than one batch here
  expect_equal(requested[1:2], c(300L, 300L))
  # every batch after round 2's first asks only for the still-missing draws,
  # never the full 300 again, and the shortfall never grows batch to batch
  round2_later <- requested[-(1:2)]
  expect_true(all(round2_later < 300L))
  expect_true(all(diff(requested[-1]) <= 0))
})

test_that("npe_sequential checks its counts before the first round", {
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + stats::rnorm(1, sd = 0.5)
  seq_call <- function(n_rounds = 2, n_simulations = 100, ...) {
    npe_sequential(prior, simulator, x_obs = 0, n_rounds = n_rounds,
                   n_simulations = n_simulations,
                   density_estimator = "linear_gaussian", ...)
  }
  expect_error(seq_call(n_rounds = 2.5),
               "`n_rounds` must be a single whole number of at least 1 since")
  expect_error(seq_call(n_simulations = 1),
               "`n_simulations` must be whole numbers of at least 2")
  expect_error(seq_call(n_simulations = c(500, 0)),
               "not 500, 0")
  expect_error(seq_call(n_truncation_samples = 0),
               "`n_truncation_samples` must be")
  expect_error(seq_call(max_proposal_batches = 0),
               "`max_proposal_batches` must be a single whole number of at least 1 since")
})

test_that("npe_sequential rejects an n_simulations vector whose length doesn't match n_rounds", {
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + stats::rnorm(1, sd = 0.5)
  # two budgets for three rounds used to be silently recycled by rep_len()
  # into c(500, 1000, 500) instead of erroring (issue #191)
  expect_error(
    npe_sequential(prior, simulator, x_obs = 0, n_rounds = 3,
                   n_simulations = c(500, 1000),
                   density_estimator = "linear_gaussian"),
    "`n_simulations` must be length 1 or 3.*not 2")
})

test_that("npe_sequential accepts n_simulations as a scalar or as a full-length vector", {
  set.seed(25)
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + stats::rnorm(1, sd = 0.5)
  scalar_fit <- npe_sequential(prior, simulator, x_obs = 0, n_rounds = 2,
                               n_simulations = 100,
                               density_estimator = "linear_gaussian", seed = 1)
  expect_equal(scalar_fit$n_simulations, 200L)
  vector_fit <- npe_sequential(prior, simulator, x_obs = 0, n_rounds = 2,
                               n_simulations = c(150, 50),
                               density_estimator = "linear_gaussian", seed = 2)
  expect_equal(vector_fit$rounds[[1]]$n_new, 150L)
  expect_equal(vector_fit$rounds[[2]]$n_new, 50L)
})
