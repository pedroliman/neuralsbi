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
  expect_lt(c2st(analytic_draws, draws, classifier = "logistic",
                 seed = 2)$accuracy, 0.6)
  if (has_torch()) expect_lt(c2st(analytic_draws, draws, seed = 2)$accuracy, 0.6)
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
  # so with only one proposal batch allowed, most rounds fall short of budget.
  # sd = 0.3 keeps round 2's acceptance count away from the 0-vs-1 edge (a
  # tighter likelihood lands exactly on that edge now that seed forwards into
  # the inner npe() call, #215, and reseeds R's base RNG each round, #214).
  prior <- prior_uniform(low = -50, high = 50)
  simulator <- function(theta) theta + rnorm(1, sd = 0.3)
  expect_warning(
    fit <- npe_sequential(prior, simulator, x_obs = 0.2, n_rounds = 2,
                          n_simulations = 500, epsilon = 0.5,
                          max_proposal_batches = 1,
                          density_estimator = "linear_gaussian", seed = 8),
    "proposal draws inside the truncated region"
  )
  expect_lt(fit$rounds[[2]]$n_new, 500L)
})

test_that("npe_sequential drops a NaN log_prob candidate instead of passing it to the simulator (#246)", {
  set.seed(50)
  # Same technique as test-posterior-nonfinite-de-draw.R: mock the generic
  # de_log_prob() -- called by log_prob.nsbi_posterior() -- to inject a NaN
  # the way an under-trained MAF/NSF would for an out-of-distribution prior
  # draw. Call 1 is the round-2 threshold's reference sample
  # (log_prob(post, ref, ...)), which stats::quantile() cannot tolerate a NaN
  # in; call 2 is the first proposal batch's log_prob(post, cand, ...), where
  # the corrupted row belongs.
  prior <- prior_normal(mean = 0, sd = 2)
  simulator <- function(theta) {
    # Before the fix, an NA-filled row from the buggy `keep` index survives
    # into theta_new and is handed to this simulator; stopping here turns
    # that into a hard test failure instead of a silently wasted simulation.
    if (anyNA(theta)) stop("simulator called with a non-finite theta")
    theta + rnorm(length(theta), sd = 0.3)
  }
  real_de_log_prob <- de_log_prob
  call_count <- 0L
  local_mocked_bindings(
    de_log_prob = function(de, theta, x) {
      call_count <<- call_count + 1L
      lp <- real_de_log_prob(de, theta, x)
      if (call_count == 2L) lp[1] <- NaN
      lp
    }
  )

  fit <- npe_sequential(prior, simulator, x_obs = 0.5, n_rounds = 2,
                        n_simulations = c(300, 30), epsilon = 0.3,
                        density_estimator = "linear_gaussian", seed = 51)
  expect_gt(call_count, 1L)
  expect_equal(fit$rounds[[2]]$n_new, 30L)
  expect_false(anyNA(fit$theta))
  # acceptance is reported against usable draws, so it cannot exceed 1 even
  # though one proposal row was rejected for a reason other than truncation
  expect_lte(fit$rounds[[2]]$acceptance, 1)
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

test_that("npe_sequential rejects an out-of-range epsilon before round 1 simulates (#226)", {
  # epsilon only gets read starting in round 2 (it sets the truncation
  # threshold), so an out-of-range value used to pass round 1 silently and
  # only fail once round 2 called stats::quantile() with the raw base-R
  # error "'probs' outside [0,1]" -- burning round 1's simulation budget on a
  # value that was never going to work.
  prior <- prior_normal(mean = 0, sd = 1)
  n_calls <- 0L
  simulator <- function(theta) {
    n_calls <<- n_calls + 1L
    theta + stats::rnorm(1, sd = 0.5)
  }
  seq_call <- function(epsilon) {
    npe_sequential(prior, simulator, x_obs = 0, n_rounds = 2,
                   n_simulations = 100, epsilon = epsilon,
                   density_estimator = "linear_gaussian")
  }
  expect_error(seq_call(1.5), "`epsilon` must be a single number")
  expect_error(seq_call(-0.1), "`epsilon` must be a single number")
  expect_error(seq_call(NA_real_), "`epsilon` must be a single number")
  expect_equal(n_calls, 0L)

  fit <- seq_call(0.1)
  expect_s3_class(fit, "nsbi_snpe")
})

test_that("npe_sequential checks round 1's estimator/training args before round 1 simulates (#251)", {
  # n_bins, batch_size and the rest of npe()'s architecture/training-control
  # arguments reach round 1 through `...` and used to be checked only once
  # the inner npe() call ran, at the end of the round -- after
  # prepare_simulations() had already spent round 1's simulation budget. The
  # same class of bug #226 fixed for n_rounds/n_simulations/epsilon/etc.
  prior <- prior_normal(mean = 0, sd = 1)
  n_calls <- 0L
  simulator <- function(theta) {
    n_calls <<- n_calls + 1L
    theta + stats::rnorm(1, sd = 0.5)
  }
  seq_call <- function(...) {
    npe_sequential(prior, simulator, x_obs = 0, n_rounds = 2,
                   n_simulations = 500, density_estimator = "linear_gaussian",
                   ...)
  }
  expect_error(seq_call(n_bins = 1), "`n_bins` must be a single whole number of at least 2")
  expect_equal(n_calls, 0L)
  expect_error(seq_call(batch_size = 0), "`batch_size` must be a single whole number of at least 1")
  expect_equal(n_calls, 0L)
  expect_error(seq_call(device = "not-a-device"), "`device` must be one of")
  expect_equal(n_calls, 0L)

  fit <- seq_call(n_bins = 4)
  expect_s3_class(fit, "nsbi_snpe")
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

test_that("npe_sequential() is reproducible given `seed`, including torch's RNG (#215)", {
  skip_if_no_torch()
  # Regression test for GitHub #215. npe_sequential() seeds R's base RNG
  # itself (set.seed(seed), above), but until the fix it never forwarded
  # `seed` to the inner npe() call (R/sequential.R), so npe()'s own `seed`
  # stayed NULL and train_restarts() (R/train.R) never called
  # torch::torch_manual_seed(): torch's RNG, which drives network weight
  # initialization every round, was never seeded by npe_sequential().
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + stats::rnorm(1, sd = 0.5)

  # Diverge torch's ambient RNG by a different amount before each call, so
  # an identical result proves it came from `seed`, not from the RNG
  # already happening to agree.
  torch::torch_manual_seed(100)
  fit1 <- npe_sequential(prior, simulator, x_obs = 0.3, n_rounds = 2,
                         n_simulations = 100, density_estimator = "mdn",
                         n_components = 1L, hidden = c(4L), max_epochs = 3L,
                         n_restarts = 1L, seed = 1)
  torch::torch_manual_seed(200)
  fit2 <- npe_sequential(prior, simulator, x_obs = 0.3, n_rounds = 2,
                         n_simulations = 100, density_estimator = "mdn",
                         n_components = 1L, hidden = c(4L), max_epochs = 3L,
                         n_restarts = 1L, seed = 1)

  expect_identical(fit1$de$history, fit2$de$history)
  expect_identical(fit1$de$best_val_loss, fit2$de$best_val_loss)
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
