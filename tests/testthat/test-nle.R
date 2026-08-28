gauss_prior <- function() {
  prior_uniform(c(mu = -3, nu = -3), c(mu = 3, nu = 3))
}

gauss_sim <- function(mu, nu) {
  c(a = mu + stats::rnorm(1, sd = 0.4), b = nu + 0.5 * mu + stats::rnorm(1, sd = 0.3))
}

lingauss_fit <- function(n = 2000, seed = 1) {
  nle(gauss_prior(), gauss_sim, n_simulations = n,
      density_estimator = "linear_gaussian", seed = seed)
}

test_that("nle() returns a fit carrying the fields the pipeline expects", {
  fit <- lingauss_fit()

  expect_s3_class(fit, "nsbi_nle")
  expect_equal(fit$dim_theta, 2L)
  expect_equal(fit$dim_x, 2L)
  expect_equal(fit$param_names, c("mu", "nu"))
  expect_equal(fit$x_names, c("a", "b"))
  expect_equal(fit$n_simulations, 2000L)
  expect_s3_class(fit$std_theta, "nsbi_standardizer")
  expect_s3_class(fit$std_x, "nsbi_standardizer")
  expect_output(print(fit), "Neural Likelihood Estimation")
  expect_output(print(fit), "q\\(x \\| theta\\)")
})

test_that("nle() learns q(x | theta), not q(theta | x)", {
  # The whole design rests on the roles being swapped relative to npe(). If
  # they were not, the estimator's target would be the parameter and the
  # density would not respond to moving x.
  fit <- lingauss_fit()
  theta <- c(1, -0.5)

  near <- log_lik(fit, theta, matrix(c(1, 0), nrow = 1))
  far <- log_lik(fit, theta, matrix(c(1, 8), nrow = 1))

  expect_gt(near, far)
})

test_that("the surrogate likelihood tracks the true one on a Gaussian model", {
  # x | theta is Gaussian and linear in theta here, so linear_gaussian is
  # exact up to estimation error and the fit can be scored against the truth.
  fit <- lingauss_fit(n = 6000, seed = 2)
  set.seed(3)
  theta <- matrix(stats::runif(20, -2, 2), ncol = 2)
  x <- matrix(c(0.3, -0.4), nrow = 1)

  Sigma <- matrix(c(0.4^2, 0, 0, 0.3^2), 2, 2)
  truth <- vapply(seq_len(nrow(theta)), function(i) {
    mu <- c(theta[i, 1], theta[i, 2] + 0.5 * theta[i, 1])
    dmvnorm_chol(x, mu, chol(Sigma))
  }, numeric(1))

  expect_equal(log_lik(fit, theta, x), truth, tolerance = 0.05)
})

test_that("independent observations add up", {
  fit <- lingauss_fit()
  theta <- rbind(c(0.5, -0.5), c(1, 1))
  x1 <- matrix(c(0.4, 0.1), nrow = 1)
  x2 <- matrix(c(-0.2, 0.7), nrow = 1)

  expect_equal(log_lik(fit, theta, rbind(x1, x2)),
               log_lik(fit, theta, x1) + log_lik(fit, theta, x2))
})

test_that("sum_iid = FALSE returns the per-observation matrix", {
  fit <- lingauss_fit()
  theta <- rbind(c(0.5, -0.5), c(1, 1), c(-1, 0))
  x <- matrix(stats::rnorm(8), ncol = 2)

  per_obs <- log_lik(fit, theta, x, sum_iid = FALSE)

  expect_equal(dim(per_obs), c(3L, 4L))
  expect_equal(rowSums(per_obs), log_lik(fit, theta, x))
})

test_that("the shape is right at every corner of the cross product", {
  # The estimators that factorize the i.i.d. sum build their result a row at a
  # time, and R drops a one-row matrix to a vector, which silently transposes
  # the answer. Every combination of one and several is checked for both.
  for (estimator in c("linear_gaussian", "mdn")) {
    if (estimator == "mdn") skip_if_no_torch()
    fit <- nle(gauss_prior(), gauss_sim, n_simulations = 600,
               density_estimator = estimator, hidden = c(16L, 16L),
               n_components = 3L, max_epochs = 10L, seed = 1)
    for (n_theta in c(1L, 5L)) {
      for (n_obs in c(1L, 3L)) {
        theta <- matrix(stats::runif(2 * n_theta, -2, 2), ncol = 2)
        x <- matrix(stats::rnorm(2 * n_obs), ncol = 2)

        per_obs <- log_lik(fit, theta, x, sum_iid = FALSE)
        label <- sprintf("%s, %d theta, %d obs", estimator, n_theta, n_obs)

        expect_equal(dim(per_obs), c(n_theta, n_obs), label = label)
        expect_equal(rowSums(per_obs), log_lik(fit, theta, x), label = label)
        # The general path is the reference: the fast ones must agree with it.
        expect_equal(
          per_obs,
          de_log_lik_iid.default(fit$de,
                                 apply_standardizer(fit$std_x, x),
                                 apply_standardizer(fit$std_theta, theta)) +
            standardizer_log_jac(fit$std_x),
          tolerance = 1e-5, ignore_attr = TRUE, label = label)
      }
    }
  }
})

test_that("blocking does not change the answer", {
  fit <- lingauss_fit()
  theta <- matrix(stats::runif(40, -2, 2), ncol = 2)
  x <- matrix(stats::rnorm(20), ncol = 2)

  expect_equal(log_lik(fit, theta, x, max_batch = 7),
               log_lik(fit, theta, x, max_batch = 1e5))
})

test_that("the blocked cross product visits every (theta, x) pair once", {
  # The check above cannot see a blocking bug: linear_gaussian has its own
  # de_log_lik_iid() and de_iid_evaluator() methods that score the whole
  # observation set per parameter and ignore max_batch entirely. cross_iid() is
  # the path every flow takes, and without torch nothing else in the suite
  # reaches it. A scorer that encodes which pair it was handed pins the layout
  # exactly, rather than letting a transposition hide inside a row sum.
  score <- function(de, x, theta) x[, 1] * 1000 + theta[, 1]
  n_obs <- 13L
  n_theta <- 37L
  x <- matrix(seq_len(n_obs), ncol = 1)
  theta <- matrix(seq_len(n_theta), ncol = 1)
  expected <- outer(seq_len(n_theta), seq_len(n_obs),
                    function(i, j) j * 1000 + i)

  # Batch sizes either side of one observation set, of two, and of the whole
  # cross product, so the short final block and the single-row block both run.
  for (max_batch in c(1, 2, 5, 13, 14, 26, 40, 1e5)) {
    label <- paste("max_batch", max_batch)
    expect_equal(iid_matrix(NULL, x, theta, max_batch, score), expected,
                 label = label)
    ev <- iid_evaluator(NULL, x, max_batch, score)
    expect_equal(ev(theta), rowSums(expected), label = label)
  }
})

test_that("blocking does not change a neural estimator's answer either", {
  # The MDN chunks over observations and the flow over parameters, and both
  # reuse a constant block that a short final block has to truncate. A batch
  # size that divides neither count exercises that truncation.
  skip_if_no_torch()
  for (estimator in c("mdn", "maf")) {
    fit <- nle(gauss_prior(), gauss_sim, n_simulations = 600,
               density_estimator = estimator, hidden = c(16L, 16L),
               n_components = 3L, n_transforms = 2L, max_epochs = 10L, seed = 1)
    theta <- matrix(stats::runif(14, -2, 2), ncol = 2)
    x <- matrix(stats::rnorm(14), ncol = 2)

    expect_equal(log_lik(fit, theta, x, max_batch = 11),
                 log_lik(fit, theta, x, max_batch = 1e5),
                 tolerance = 1e-6, label = estimator)
    expect_equal(log_lik(fit, theta, x, sum_iid = FALSE, max_batch = 11),
                 log_lik(fit, theta, x, sum_iid = FALSE, max_batch = 1e5),
                 tolerance = 1e-6, label = estimator)
  }
})

test_that("mdn_theta_chunk_size() bounds theta the same way cross_iid() bounds it (#240)", {
  # Pure arithmetic, no torch needed: theta_chunk * n_obs stays within
  # max_batch, and never exceeds n_theta itself.
  expect_equal(mdn_theta_chunk_size(n_theta = 1000, n_obs = 1, max_batch = 50), 50L)
  expect_equal(mdn_theta_chunk_size(n_theta = 1000, n_obs = 200, max_batch = 50), 1L)
  expect_equal(mdn_theta_chunk_size(n_theta = 10, n_obs = 1, max_batch = 1e5), 10L)
})

test_that("the MDN i.i.d. path chunks theta, not just x, under max_batch (#240)", {
  # Before the fix, mdn_iid_blocks() ran mdn_mixture() -- the MLP forward pass
  # and Cholesky assembly -- over the whole of theta in one call regardless of
  # max_batch, only chunking the observation side. A dense theta grid against a
  # handful of observations (e.g. a profile-likelihood plot) is exactly the
  # shape that bypassed max_batch entirely: n_theta large, n_obs small.
  skip_if_no_torch()
  fit <- nle(gauss_prior(), gauss_sim, n_simulations = 600,
             density_estimator = "mdn", hidden = c(16L, 16L),
             n_components = 3L, max_epochs = 10L, seed = 1)
  theta <- matrix(stats::runif(2 * 47, -2, 2), ncol = 2)  # 47 rows of theta
  x <- matrix(stats::rnorm(4), ncol = 2)                  # 2 observations

  orig_mixture <- mdn_mixture
  sizes <- integer(0)
  local_mocked_bindings(mdn_mixture = function(de, tt) {
    sizes <<- c(sizes, tt$shape[1])
    orig_mixture(de, tt)
  })

  chunked <- log_lik(fit, theta, x, sum_iid = FALSE, max_batch = 10)
  # theta chunk size for n_obs = 2, max_batch = 10 is floor(10 / 2) = 5, so 47
  # rows of theta take 10 calls to mdn_mixture(), the last one short -- not the
  # single call the bug produced regardless of max_batch.
  expect_equal(length(sizes), 10L)
  expect_true(all(sizes <= 5L))
  expect_equal(sum(sizes), 47L)

  unchunked <- log_lik(fit, theta, x, sum_iid = FALSE, max_batch = 1e5)
  expect_equal(chunked, unchunked, tolerance = 1e-6)

  # The summed evaluator (log_lik(sum_iid = TRUE), and every MCMC step) goes
  # through the same mdn_iid_blocks(), so it has to agree too.
  expect_equal(log_lik(fit, theta, x, max_batch = 10),
               log_lik(fit, theta, x, max_batch = 1e5), tolerance = 1e-6)
})

test_that("the traced MDN likelihood agrees with the eager one", {
  # de_iid_evaluator() serves the first few calls eagerly and then records a
  # TorchScript trace per parameter-row count. The trace is checked internally
  # before use, but that check is the thing most worth checking, so this walks
  # past the warmup at two different row counts and compares every call.
  skip_if_no_torch()
  # withr is a Suggests, so under _R_CHECK_FORCE_SUGGESTS_=false this has to
  # skip rather than error on the with_options() calls below.
  skip_if_not_installed("withr")
  fit <- nle(gauss_prior(), gauss_sim, n_simulations = 800,
             density_estimator = "mdn", hidden = c(16L, 16L),
             n_components = 3L, max_epochs = 10L, seed = 1)
  x_z <- apply_standardizer(fit$std_x, matrix(stats::rnorm(30), ncol = 2))
  theta_z <- apply_standardizer(fit$std_theta,
                                matrix(stats::runif(20, -2, 2), ncol = 2))
  one <- theta_z[3, , drop = FALSE]

  ref_many <- withr::with_options(
    list(neuralsbi.jit = FALSE), de_iid_evaluator(fit$de, x_z)(theta_z))
  ref_one <- withr::with_options(
    list(neuralsbi.jit = FALSE), de_iid_evaluator(fit$de, x_z)(one))

  ev <- de_iid_evaluator(fit$de, x_z)
  for (i in seq_len(8)) {
    expect_equal(ev(theta_z), ref_many, tolerance = 1e-6)
    expect_equal(ev(one), ref_one, tolerance = 1e-6)
  }
})

test_that("the summed path agrees with summing the per-observation matrix", {
  # log_lik(sum_iid = TRUE) does not build the n_theta x n_obs matrix and sum
  # it; it goes through de_iid_evaluator(), which sums where the densities are
  # produced. The two have to agree, including at the corners where one of the
  # counts is 1 and a matrix would collapse to a vector.
  for (estimator in c("linear_gaussian", "mdn", "maf")) {
    if (estimator != "linear_gaussian") skip_if_no_torch()
    fit <- nle(gauss_prior(), gauss_sim, n_simulations = 600,
               density_estimator = estimator, hidden = c(16L, 16L),
               n_components = 3L, n_transforms = 2L, max_epochs = 10L, seed = 1)
    for (n_theta in c(1L, 4L)) {
      for (n_obs in c(1L, 5L)) {
        theta <- matrix(stats::runif(2 * n_theta, -2, 2), ncol = 2)
        x <- matrix(stats::rnorm(2 * n_obs), ncol = 2)
        label <- sprintf("%s, %d theta, %d obs", estimator, n_theta, n_obs)

        expect_equal(log_lik(fit, theta, x),
                     rowSums(log_lik(fit, theta, x, sum_iid = FALSE)),
                     tolerance = 1e-6, label = label)
      }
    }
  }
})

test_that("wrongly shaped input is rejected by name", {
  fit <- lingauss_fit()

  expect_error(log_lik(fit, matrix(0, 2, 3), matrix(0, 1, 2)),
               "`theta` must have 2 columns")
  expect_error(log_lik(fit, matrix(0, 2, 2), matrix(0, 1, 3)),
               "`x` must have 2 columns")
})

test_that("log_lik() rejects a non-finite theta or x instead of returning NA (#202)", {
  fit <- lingauss_fit()
  theta <- matrix(c(1, -0.5), nrow = 1)
  x <- matrix(c(0.4, -0.2), nrow = 1)

  expect_error(log_lik(fit, theta, matrix(c(0.4, NA), nrow = 1)),
               "`x` contains 1 non-finite value \\(NA\\), first at row 1, column 2")
  expect_error(log_lik(fit, theta, matrix(c(NaN, -0.2), nrow = 1)),
               "`x` contains")
  expect_error(log_lik(fit, theta, matrix(c(Inf, -0.2), nrow = 1)),
               "`x` contains 1 non-finite value \\(Inf\\)")
  expect_error(log_lik(fit, matrix(c(1, NA), nrow = 1), x),
               "`theta` contains 1 non-finite value \\(NA\\), first at row 1, column 2")
})

test_that("log_lik() rejects a non-finite max_batch instead of an opaque rep() error (#230)", {
  fit <- lingauss_fit()
  theta <- matrix(c(1, -0.5), nrow = 1)
  x <- matrix(c(0.4, -0.2), nrow = 1)

  expect_error(log_lik(fit, theta, x, max_batch = NA), "`max_batch`")
  expect_error(log_lik(fit, theta, x, max_batch = NaN), "`max_batch`")
})

test_that("nle() needs either a simulator or pre-computed simulations", {
  expect_error(nle(gauss_prior()), "Provide either")
})

test_that("a bad density_estimator errors before the simulator runs", {
  calls <- 0L
  counting_simulator <- function(mu, nu) {
    calls <<- calls + 1L
    gauss_sim(mu, nu)
  }
  expect_error(
    nle(gauss_prior(), counting_simulator, n_simulations = 100,
        density_estimator = "mfa"),
    "should be one of"
  )
  expect_identical(calls, 0L)
})

test_that("pre-computed simulations give the same fit as running the simulator", {
  sims <- simulate_for_sbi(gauss_sim, gauss_prior(), n = 800, seed = 9)
  from_sims <- nle(gauss_prior(), theta = sims$theta, x = sims$x,
                   density_estimator = "linear_gaussian")
  direct <- nle(gauss_prior(), gauss_sim, n_simulations = 800,
                density_estimator = "linear_gaussian", seed = 9)

  expect_equal(from_sims$de$B, direct$de$B)
})

test_that("likelihood_fn() is a vectorized closure over theta", {
  fit <- lingauss_fit()
  x_obs <- matrix(stats::rnorm(20), ncol = 2)
  loglik <- likelihood_fn(fit, x_obs)
  theta <- matrix(stats::runif(10, -2, 2), ncol = 2)

  expect_type(loglik, "closure")
  expect_equal(loglik(theta), log_lik(fit, theta, x_obs))
  expect_length(loglik(c(0, 0)), 1L)
  expect_equal(attr(loglik, "x_obs"), x_obs)
})

test_that("likelihood_fn() drops into optim() without adaptation", {
  # The point of the closure: nothing downstream needs to know about neuralsbi.
  set.seed(10)
  prior <- prior_uniform(c(mu = -4), c(mu = 4))
  fit <- nle(prior, function(mu) c(y = stats::rnorm(1, mu, 0.5)),
             n_simulations = 3000, density_estimator = "linear_gaussian",
             seed = 4)
  x_obs <- matrix(stats::rnorm(200, mean = 1.5, sd = 0.5), ncol = 1)

  mle <- stats::optimize(likelihood_fn(fit, x_obs), c(-4, 4), maximum = TRUE)

  expect_equal(mle$maximum, mean(x_obs), tolerance = 0.1)
})

test_that("a fit survives the save/load round trip", {
  fit <- lingauss_fit()
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  save_nle(fit, path)
  again <- load_nle(path)

  expect_s3_class(again, "nsbi_nle")
  expect_equal(log_lik(again, c(0.5, 0.5), matrix(c(0, 0), nrow = 1)),
               log_lik(fit, c(0.5, 0.5), matrix(c(0, 0), nrow = 1)))
})

test_that("the neural estimators fit and evaluate", {
  skip_if_no_torch()
  for (estimator in c("mdn", "maf", "nsf")) {
    fit <- nle(gauss_prior(), gauss_sim, n_simulations = 600,
               density_estimator = estimator, hidden = c(16L, 16L),
               n_components = 3L, n_transforms = 2L, max_epochs = 15L, seed = 1)
    x <- matrix(stats::rnorm(6), ncol = 2)
    lp <- log_lik(fit, rbind(c(0, 0), c(1, 1)), x)

    expect_equal(fit$density_estimator, estimator)
    expect_length(lp, 2L)
    expect_true(all(is.finite(lp)), label = estimator)
  }
})

test_that("a neural fit reproduces a Gaussian likelihood", {
  skip_if_no_torch()
  fit <- nle(gauss_prior(), gauss_sim, n_simulations = 6000,
             density_estimator = "maf", max_epochs = 300L, seed = 5)
  set.seed(6)
  theta <- matrix(stats::runif(20, -1.5, 1.5), ncol = 2)
  x <- matrix(c(0.3, -0.4), nrow = 1)

  Sigma <- matrix(c(0.4^2, 0, 0, 0.3^2), 2, 2)
  truth <- vapply(seq_len(nrow(theta)), function(i) {
    mu <- c(theta[i, 1], theta[i, 2] + 0.5 * theta[i, 1])
    dmvnorm_chol(x, mu, chol(Sigma))
  }, numeric(1))

  # A flow is not exact the way the closed-form baseline is; this checks it
  # gets the shape right, not that it nails the density.
  expect_gt(stats::cor(log_lik(fit, theta, x), truth), 0.97)
})
