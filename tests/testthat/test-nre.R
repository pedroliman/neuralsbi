gauss_prior <- function() {
  prior_uniform(c(mu = -3, nu = -3), c(mu = 3, nu = 3))
}

gauss_sim <- function(mu, nu) {
  c(a = mu + stats::rnorm(1, sd = 0.4), b = nu + 0.5 * mu + stats::rnorm(1, sd = 0.3))
}

logistic_fit <- function(n = 4000, seed = 1) {
  nre(gauss_prior(), gauss_sim, n_simulations = n, classifier = "logistic",
      seed = seed)
}

# The log ratio differs from the log likelihood by log p(x), a function of x
# alone, so a comparison against the truth is only ever up to an offset. Every
# accuracy check below removes it the same way.
centred <- function(v) v - mean(v)

test_that("nre() returns a fit carrying the fields the pipeline expects", {
  fit <- logistic_fit()

  expect_s3_class(fit, "nsbi_nre")
  expect_equal(fit$dim_theta, 2L)
  expect_equal(fit$dim_x, 2L)
  expect_equal(fit$param_names, c("mu", "nu"))
  expect_equal(fit$x_names, c("a", "b"))
  expect_equal(fit$n_simulations, 4000L)
  expect_equal(fit$classifier, "logistic")
  expect_s3_class(fit$std_theta, "nsbi_standardizer")
  expect_s3_class(fit$std_x, "nsbi_standardizer")
  expect_output(print(fit), "Neural Ratio Estimation")
  expect_output(print(fit), "p\\(x \\| theta\\) / p\\(x\\)")
  expect_output(print(fit), "contrastive atoms : 10")
})

test_that("the learned ratio tracks the true likelihood up to log p(x)", {
  # x | theta is Gaussian and linear in theta, so the quadratic feature basis
  # spans the log ratio's parameter dependence exactly and the fit can be
  # scored against the truth.
  fit <- logistic_fit(n = 8000, seed = 2)
  set.seed(3)
  theta <- matrix(stats::runif(40, -1.5, 1.5), ncol = 2)
  x <- matrix(c(0.3, -0.4), nrow = 1)

  Sigma <- matrix(c(0.4^2, 0, 0, 0.3^2), 2, 2)
  truth <- vapply(seq_len(nrow(theta)), function(i) {
    mu <- c(theta[i, 1], theta[i, 2] + 0.5 * theta[i, 1])
    dmvnorm_chol(x, mu, chol(Sigma))
  }, numeric(1))

  expect_equal(centred(log_ratio(fit, theta, x)), centred(truth),
               tolerance = 0.15)
})

test_that("the ratio responds to x, not just to theta", {
  fit <- logistic_fit()
  theta <- c(1, -0.5)

  near <- log_ratio(fit, theta, matrix(c(1, 0), nrow = 1))
  far <- log_ratio(fit, theta, matrix(c(1, 8), nrow = 1))

  expect_gt(near, far)
})

test_that("independent observations add up", {
  fit <- logistic_fit()
  theta <- rbind(c(0.5, -0.5), c(1, 1))
  x1 <- matrix(c(0.4, 0.1), nrow = 1)
  x2 <- matrix(c(-0.2, 0.7), nrow = 1)

  expect_equal(log_ratio(fit, theta, rbind(x1, x2)),
               log_ratio(fit, theta, x1) + log_ratio(fit, theta, x2))
})

test_that("sum_iid = FALSE returns the per-observation matrix", {
  fit <- logistic_fit()
  theta <- rbind(c(0.5, -0.5), c(1, 1), c(-1, 0))
  x <- matrix(stats::rnorm(8), ncol = 2)

  per_obs <- log_ratio(fit, theta, x, sum_iid = FALSE)

  expect_equal(dim(per_obs), c(3L, 4L))
  expect_equal(rowSums(per_obs), log_ratio(fit, theta, x))
})

test_that("the shape is right at every corner of the cross product", {
  fit <- logistic_fit(n = 1000)
  for (n_theta in c(1L, 5L)) {
    for (n_obs in c(1L, 3L)) {
      theta <- matrix(stats::runif(2 * n_theta, -2, 2), ncol = 2)
      x <- matrix(stats::rnorm(2 * n_obs), ncol = 2)
      label <- sprintf("%d theta, %d obs", n_theta, n_obs)

      per_obs <- log_ratio(fit, theta, x, sum_iid = FALSE)
      expect_equal(dim(per_obs), c(n_theta, n_obs), label = label)
      expect_equal(rowSums(per_obs), log_ratio(fit, theta, x), label = label)
    }
  }
})

test_that("blocking does not change the answer", {
  fit <- logistic_fit(n = 1000)
  theta <- matrix(stats::runif(40, -2, 2), ncol = 2)
  x <- matrix(stats::rnorm(20), ncol = 2)

  expect_equal(log_ratio(fit, theta, x, max_batch = 7),
               log_ratio(fit, theta, x, max_batch = 1e5))
  expect_equal(log_ratio(fit, theta, x, sum_iid = FALSE, max_batch = 7),
               log_ratio(fit, theta, x, sum_iid = FALSE, max_batch = 1e5))
})

test_that("standardizing x leaves the ratio unchanged", {
  # A ratio needs no change-of-variables term: the Jacobian cancels between
  # p(x | theta) and p(x). Rescaling the data by a constant therefore has to
  # leave the fitted log ratio alone, which the log-likelihood of an nle() fit
  # would not.
  set.seed(11)
  sims <- simulate_for_sbi(gauss_sim, gauss_prior(), n = 3000)
  plain <- nre(gauss_prior(), theta = sims$theta, x = sims$x,
               classifier = "logistic")
  scaled <- nre(gauss_prior(), theta = sims$theta, x = sims$x * 100,
                classifier = "logistic")
  theta <- matrix(stats::runif(10, -2, 2), ncol = 2)
  x <- matrix(stats::rnorm(4), ncol = 2)

  expect_equal(log_ratio(plain, theta, x),
               log_ratio(scaled, theta, x * 100), tolerance = 1e-4)
})

test_that("the posterior recovers a conjugate Gaussian", {
  # theta ~ U(-3, 3)^2 wide enough to act flat, x | theta ~ N(theta, 0.5^2),
  # 40 independent observations. The posterior is then Gaussian with mean
  # colMeans(x_obs) and sd 0.5 / sqrt(40) in each coordinate.
  set.seed(20)
  prior <- prior_uniform(c(mu = -3, nu = -3), c(mu = 3, nu = 3))
  sim <- function(mu, nu) c(a = stats::rnorm(1, mu, 0.5),
                            b = stats::rnorm(1, nu, 0.5))
  fit <- nre(prior, sim, n_simulations = 8000, classifier = "logistic",
             seed = 21)
  x_obs <- cbind(stats::rnorm(40, 0.8, 0.5), stats::rnorm(40, -0.6, 0.5))

  post <- posterior(fit, x_obs, n_chains = 8, warmup = 200, seed = 22)
  draws <- sample(post, 2000)

  expect_equal(colMeans(draws), colMeans(x_obs), tolerance = 0.05,
               ignore_attr = TRUE)
  expect_equal(unname(apply(draws, 2, stats::sd)), rep(0.5 / sqrt(40), 2),
               tolerance = 0.25)
})

test_that("the posterior object reports itself and caches its draws", {
  fit <- logistic_fit(n = 2000)
  x_obs <- matrix(stats::rnorm(20), ncol = 2)
  post <- posterior(fit, x_obs, n_chains = 4, warmup = 50, seed = 5)

  expect_s3_class(post, "nsbi_nre_posterior")
  expect_s3_class(post, "nsbi_posterior")
  expect_output(print(post), "nsbi_nre_posterior")
  # No Stan export for a ratio estimator, so the print method must not offer one.
  expect_output(print(post), "map_estimate\\(post\\)")
  expect_false(any(grepl("stan_code", utils::capture.output(print(post)))))

  first <- sample(post, 200)
  expect_equal(dim(first), c(200L, 2L))
  expect_equal(colnames(first), c("mu", "nu"))
  expect_false(is.null(attr(first, "diagnostics")))
  # Second call is served from the cache, so it is identical rather than merely
  # distributed the same way.
  expect_identical(sample(post, 200), first)
})

test_that("log_prob() on an NRE posterior is the unnormalized potential", {
  fit <- logistic_fit(n = 2000)
  x_obs <- matrix(stats::rnorm(10), ncol = 2)
  post <- posterior(fit, x_obs, n_chains = 4, warmup = 20)
  theta <- rbind(c(0.5, -0.5), c(1, 1))

  lp <- log_prob(post, theta)
  expect_equal(lp, log_ratio(fit, theta, x_obs) +
                     as.numeric(fit$prior$log_prob(theta)))
  expect_warning(log_prob(post, theta, normalize = TRUE), "no normalizing")
  expect_equal(log_prob(post, rbind(c(9, 9))), -Inf)
})

# GitHub #163: mcmc_log_prob() is shared by NLE and NRE posteriors, so the
# same fix that gives log_prob() a named error for a non-numeric theta on an
# NLE fit applies here too.
test_that("log_prob() rejects non-numeric theta by name instead of returning NA", {
  fit <- logistic_fit(n = 2000)
  x_obs <- matrix(stats::rnorm(10), ncol = 2)
  post <- posterior(fit, x_obs, n_chains = 4, warmup = 20)

  expect_error(
    log_prob(post, data.frame(mu = c("x", "y"), nu = c(1, 2)),
             normalize = FALSE),
    "`theta` has non-numeric columns: mu")
})

test_that("wrongly shaped input is rejected by name", {
  fit <- logistic_fit(n = 500)

  expect_error(log_ratio(fit, matrix(0, 2, 3), matrix(0, 1, 2)),
               "`theta` must have 2 columns")
  expect_error(log_ratio(fit, matrix(0, 2, 2), matrix(0, 1, 3)),
               "`x` must have 2 columns")
})

test_that("log_ratio() rejects a non-finite theta or x instead of returning NA (#202)", {
  fit <- logistic_fit(n = 500)
  theta <- matrix(c(1, -0.5), nrow = 1)
  x <- matrix(c(0.4, -0.2), nrow = 1)

  expect_error(log_ratio(fit, theta, matrix(c(0.4, NA), nrow = 1)),
               "`x` contains 1 non-finite value \\(NA\\), first at row 1, column 2")
  expect_error(log_ratio(fit, theta, matrix(c(NaN, -0.2), nrow = 1)),
               "`x` contains")
  expect_error(log_ratio(fit, theta, matrix(c(Inf, -0.2), nrow = 1)),
               "`x` contains 1 non-finite value \\(Inf\\)")
  expect_error(log_ratio(fit, matrix(c(1, NA), nrow = 1), x),
               "`theta` contains 1 non-finite value \\(NA\\), first at row 1, column 2")
})

test_that("log_ratio() rejects a non-finite max_batch instead of an opaque rep() error (#230)", {
  fit <- logistic_fit(n = 500)
  theta <- matrix(c(1, -0.5), nrow = 1)
  x <- matrix(c(0.4, -0.2), nrow = 1)

  expect_error(log_ratio(fit, theta, x, max_batch = NA), "`max_batch`")
  expect_error(log_ratio(fit, theta, x, max_batch = NaN), "`max_batch`")
})

test_that("nre() checks its arguments before the simulator runs", {
  calls <- 0L
  counting_simulator <- function(mu, nu) {
    calls <<- calls + 1L
    gauss_sim(mu, nu)
  }
  expect_error(
    nre(gauss_prior(), counting_simulator, n_simulations = 100,
        classifier = "resnest"),
    "should be one of")
  expect_error(
    nre(gauss_prior(), counting_simulator, n_simulations = 100,
        classifier = "logistic", num_atoms = 1L),
    "`num_atoms` must be a single whole number of at least 2")
  expect_error(
    nre(gauss_prior(), counting_simulator, n_simulations = 100,
        classifier = "logistic", hidden = 0L),
    "`hidden` must be a single whole number of at least 1")
  expect_error(
    nre(gauss_prior(), counting_simulator, n_simulations = 100,
        classifier = "logistic", n_blocks = -1L),
    "`n_blocks` must be a single whole number")
  expect_identical(calls, 0L)
})

test_that("nre() needs either a simulator or pre-computed simulations", {
  # logistic, so torch availability (check_torch_for_estimator(), #250) has
  # nothing to say and this is purely about prepare_simulations().
  expect_error(nre(gauss_prior(), classifier = "logistic"), "Provide either")
})

# GitHub #188: a validation split of exactly one row made the atomic
# contrastive objective return a constant zero loss every epoch, which broke
# early stopping silently -- the fit stopped after `patience` epochs holding
# the untrained epoch-1 network, reporting best_val_loss = 0 as if it were a
# perfect fit. These are pure argument-validation checks: no torch involved,
# since check_train_controls() (which train_conditional_de() calls before
# ever touching torch) is what actually enforces the floor.
test_that("nre() fails before simulating rather than train on a 1-row validation split", {
  calls <- 0L
  counting_simulator <- function(mu, nu) {
    calls <<- calls + 1L
    gauss_sim(mu, nu)
  }
  # validation_fraction defaults to 0.1, so n_simulations = 15 leaves exactly
  # 1 validation row -- the trigger reported in #188.
  expect_error(
    nre(gauss_prior(), counting_simulator, n_simulations = 15,
        classifier = "resnet"),
    "needs at least .* to score its objective on")
  expect_identical(calls, 0L)

  # The closed-form logistic classifier never splits off a validation set, so
  # the same n_simulations must not trip this floor for it.
  expect_no_error(
    nre(gauss_prior(), counting_simulator, n_simulations = 15,
        classifier = "logistic"))
  expect_identical(calls, 15L)
})

test_that("fit_nre_net() rejects a validation split too small for the atomic loss", {
  theta <- matrix(stats::rnorm(15), ncol = 1)
  x <- matrix(stats::rnorm(15), ncol = 1)

  expect_error(fit_nre_net(theta, x, classifier = "resnet"),
               "needs at least .* to score its objective on")
})

# GitHub #239: check_train_controls() enforced min_val_rows on the validation
# side only. A large validation_fraction can hold out enough rows to clear
# that floor while leaving fewer than min_val_rows for training -- e.g.
# n_simulations = 4, validation_fraction = 0.75 gives n_val = 3 (clears the
# floor of 2) and n_tr = 1 (does not). Before the fix this slipped through:
# every training step scored nre_atomic_log_prob()'s k < 2L branch, a
# constant zero loss with no gradient, and training ran silently to
# `patience` epochs reporting a best_val_loss as if it had actually trained.
test_that("nre() fails before simulating rather than train on a 1-row training split", {
  calls <- 0L
  counting_simulator <- function(mu, nu) {
    calls <<- calls + 1L
    gauss_sim(mu, nu)
  }
  expect_error(
    nre(gauss_prior(), counting_simulator, n_simulations = 4,
        validation_fraction = 0.75, classifier = "resnet"),
    "needs at least .* to score its objective on")
  expect_identical(calls, 0L)

  # The closed-form logistic classifier never splits off a validation set
  # (min_val_rows = 1L), so the same split must not trip this floor for it.
  expect_no_error(
    nre(gauss_prior(), counting_simulator, n_simulations = 4,
        validation_fraction = 0.75, classifier = "logistic"))
  expect_identical(calls, 4L)
})

test_that("fit_nre_net() rejects a training split too small for the atomic loss", {
  theta <- matrix(stats::rnorm(4), ncol = 1)
  x <- matrix(stats::rnorm(4), ncol = 1)

  expect_error(
    fit_nre_net(theta, x, classifier = "resnet", validation_fraction = 0.75),
    "needs at least .* to score its objective on")
})

test_that("nre() defers to prepare_simulations() when n_simulations can't hint a row count", {
  # An invalid n_simulations (not >= 1) and no pre-computed theta/x means the
  # early min_val_rows check has nothing to check against yet, so it must not
  # error itself -- the existing n_simulations validation inside
  # prepare_simulations() is still what reports this. classifier = "logistic"
  # keeps torch availability (check_torch_for_estimator(), #250) out of it,
  # since that is an independent check this test is not about.
  expect_error(
    nre(gauss_prior(), gauss_sim, n_simulations = 0, classifier = "logistic"),
    "`n_simulations`")
})

test_that("an embedding net is rejected by the logistic classifier", {
  expect_warning(
    nre(gauss_prior(), gauss_sim, n_simulations = 200, classifier = "logistic",
        embedding_net = embedding_mlp(output_dim = 2L)),
    "ignored by the logistic classifier")
  expect_error(
    nre(gauss_prior(), gauss_sim, n_simulations = 200, embedding_net = "no"),
    "must be built with embedding_mlp")
})

test_that("the logistic fit is deterministic given the same simulations", {
  # Contrasts are cyclic shifts, not random draws, so no RNG enters the fit.
  sims <- simulate_for_sbi(gauss_sim, gauss_prior(), n = 800, seed = 9)
  a <- nre(gauss_prior(), theta = sims$theta, x = sims$x,
           classifier = "logistic")
  b <- nre(gauss_prior(), theta = sims$theta, x = sims$x,
           classifier = "logistic")

  expect_identical(a$de$beta, b$de$beta)
})

test_that("a custom classifier is accepted as a function", {
  sims <- simulate_for_sbi(gauss_sim, gauss_prior(), n = 500, seed = 8)
  fit <- nre(gauss_prior(), theta = sims$theta, x = sims$x,
             classifier = function(theta, x) fit_logistic_ratio(theta, x))

  expect_equal(fit$classifier, "custom")
  expect_length(log_ratio(fit, rbind(c(0, 0), c(1, 1)),
                          matrix(c(0, 0), nrow = 1)), 2L)
})

test_that("a fit survives the save/load round trip", {
  fit <- logistic_fit(n = 1000)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  save_nre(fit, path)
  again <- load_nre(path)

  expect_s3_class(again, "nsbi_nre")
  expect_equal(log_ratio(again, c(0.5, 0.5), matrix(c(0, 0), nrow = 1)),
               log_ratio(fit, c(0.5, 0.5), matrix(c(0, 0), nrow = 1)))
})

test_that("summary() reports the classifier rather than a density estimator", {
  fit <- logistic_fit(n = 500)
  info <- expect_output(summary(fit), "Neural Ratio Estimation")

  expect_equal(info$classifier, "logistic")
  expect_equal(info$dim_theta, 2L)
  expect_null(info$density_estimator)
})

test_that("sbc() accepts an NRE fit", {
  set.seed(30)
  fit <- logistic_fit(n = 2000, seed = 31)
  res <- sbc(fit, gauss_sim, n_sbc = 4L, n_posterior_samples = 60L,
             n_chains = 4L, warmup = 30L)

  expect_s3_class(res, "nsbi_sbc")
  expect_equal(dim(res$ranks), c(4L, 2L))
})

test_that("the atomic row index picks the true parameter and distinct contrasts", {
  set.seed(4)
  b <- 12L
  k <- 4L
  rows <- matrix(neuralsbi:::nre_atom_rows(b, k), nrow = k)

  expect_equal(rows[1, ], seq_len(b))
  for (i in seq_len(b)) {
    expect_length(unique(rows[, i]), k)
    expect_false(i %in% rows[-1L, i])
  }
})

test_that("the atomic row index can be frozen without disturbing the stream", {
  set.seed(7)
  before <- stats::runif(1)
  set.seed(7)
  invisible(stats::runif(1))
  a <- neuralsbi:::nre_atom_rows(6L, 3L, deterministic = TRUE)
  after_frozen <- stats::runif(1)

  set.seed(7)
  invisible(stats::runif(1))
  b <- neuralsbi:::nre_atom_rows(6L, 3L, deterministic = TRUE)
  after_again <- stats::runif(1)

  expect_identical(a, b)
  expect_identical(after_frozen, after_again)
  expect_false(identical(before, after_frozen))
})

test_that("the neural classifiers fit and evaluate", {
  skip_if_no_torch()
  for (kind in c("resnet", "mlp", "linear")) {
    fit <- nre(gauss_prior(), gauss_sim, n_simulations = 600,
               classifier = kind, hidden = 16L, n_blocks = 2L,
               max_epochs = 15L, seed = 1)
    x <- matrix(stats::rnorm(6), ncol = 2)
    lr <- log_ratio(fit, rbind(c(0, 0), c(1, 1)), x)

    expect_equal(fit$classifier, kind)
    expect_length(lr, 2L)
    expect_true(all(is.finite(lr)), label = kind)
    expect_true(is.finite(fit$de$best_val_loss), label = kind)
  }
})

test_that("a neural fit survives the save/load round trip", {
  skip_if_no_torch()
  fit <- nre(gauss_prior(), gauss_sim, n_simulations = 600,
             classifier = "resnet", hidden = 16L, max_epochs = 10L, seed = 1)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  save_nre(fit, path)
  again <- load_nre(path)

  expect_equal(log_ratio(again, rbind(c(0.5, 0.5), c(-1, 0.2)),
                         matrix(c(0, 0), nrow = 1)),
               log_ratio(fit, rbind(c(0.5, 0.5), c(-1, 0.2)),
                         matrix(c(0, 0), nrow = 1)))
})

test_that("blocking does not change a neural classifier's answer either", {
  skip_if_no_torch()
  fit <- nre(gauss_prior(), gauss_sim, n_simulations = 600, hidden = 16L,
             max_epochs = 10L, seed = 1)
  theta <- matrix(stats::runif(14, -2, 2), ncol = 2)
  x <- matrix(stats::rnorm(14), ncol = 2)

  expect_equal(log_ratio(fit, theta, x, max_batch = 11),
               log_ratio(fit, theta, x, max_batch = 1e5), tolerance = 1e-6)
  expect_equal(log_ratio(fit, theta, x, sum_iid = FALSE, max_batch = 11),
               log_ratio(fit, theta, x, sum_iid = FALSE, max_batch = 1e5),
               tolerance = 1e-6)
})

test_that("a neural classifier reproduces a Gaussian likelihood ratio", {
  skip_if_no_torch()
  fit <- nre(gauss_prior(), gauss_sim, n_simulations = 6000,
             max_epochs = 300L, seed = 5)
  set.seed(6)
  theta <- matrix(stats::runif(40, -1.5, 1.5), ncol = 2)
  x <- matrix(c(0.3, -0.4), nrow = 1)

  Sigma <- matrix(c(0.4^2, 0, 0, 0.3^2), 2, 2)
  truth <- vapply(seq_len(nrow(theta)), function(i) {
    mu <- c(theta[i, 1], theta[i, 2] + 0.5 * theta[i, 1])
    dmvnorm_chol(x, mu, chol(Sigma))
  }, numeric(1))

  # A trained classifier is not exact the way the closed-form baseline is; this
  # checks it gets the shape right, not that it nails the ratio. The bar is
  # looser than the MAF's 0.97 in test-nle.R on purpose: at this budget the
  # correlation runs about 0.96 to 0.98 across seeds, and climbs with more
  # simulations (0.978 at 20k), so a 0.97 bar would fail on some runs for no
  # reason but the draw.
  expect_gt(stats::cor(log_ratio(fit, theta, x), truth), 0.95)
})

test_that("a neural classifier's posterior samples", {
  # The classifiers are checked above and the MCMC path below the logistic
  # baseline; this is the one test that runs both together.
  skip_if_no_torch()
  set.seed(40)
  prior <- prior_uniform(c(mu = -3), c(mu = 3))
  fit <- nre(prior, function(mu) c(y = stats::rnorm(1, mu, 0.5)),
             n_simulations = 4000, hidden = 32L, max_epochs = 100L, seed = 41)
  x_obs <- matrix(stats::rnorm(20, 0.7, 0.5), ncol = 1)

  post <- posterior(fit, x_obs, n_chains = 8, warmup = 100, seed = 42)
  draws <- sample(post, 1000)

  expect_equal(dim(draws), c(1000L, 1L))
  expect_true(all(is.finite(draws)))
  expect_equal(mean(draws), mean(x_obs), tolerance = 0.2)
})

test_that("an embedding net trains jointly with the classifier", {
  skip_if_no_torch()
  fit <- nre(gauss_prior(), gauss_sim, n_simulations = 600, hidden = 16L,
             max_epochs = 10L, seed = 1,
             embedding_net = embedding_mlp(output_dim = 4L, hidden = 8L))

  expect_equal(fit$de$embedding$output_dim, 4L)
  expect_output(print(fit), "embedding \\(mlp\\)")
  expect_true(all(is.finite(log_ratio(fit, rbind(c(0, 0), c(1, 1)),
                                      matrix(c(0.2, -0.1), nrow = 1)))))
})

test_that("the logistic ratio does not care what scale the data arrive on", {
  # The ridge is relative to each feature's own scale, so the fit is the same
  # whether or not the pipeline standardized first. Without that, a simulator
  # whose output has sd 5e-4 fits an absolute 1e-6 ridge against quadratic
  # features of order 1e-13 and the coefficients collapse to nothing.
  prior <- prior_uniform(c(mu = -3e-3), c(mu = 3e-3))
  sims <- simulate_for_sbi(function(mu) c(y = stats::rnorm(1, mu, 5e-4)),
                           prior, n = 3000, seed = 3)
  theta <- matrix(seq(-2e-3, 2e-3, length.out = 9), ncol = 1)
  x <- matrix(1e-3, nrow = 1)
  truth <- stats::dnorm(1e-3, theta[, 1], 5e-4, log = TRUE)

  raw <- nre(prior, theta = sims$theta, x = sims$x, classifier = "logistic",
             standardize = FALSE)
  scaled <- nre(prior, theta = sims$theta, x = sims$x, classifier = "logistic",
                standardize = TRUE)

  expect_equal(centred(log_ratio(raw, theta, x)), centred(truth),
               tolerance = 0.05)
  # Up to the offset, and only up to it: the features that depend on `x` alone
  # difference away during the fit, so their coefficients are whatever the
  # ridge leaves them and they do not survive a change of units. That is the
  # part of the log ratio the atomic objective never pins down anyway.
  expect_equal(centred(log_ratio(raw, theta, x)),
               centred(log_ratio(scaled, theta, x)), tolerance = 1e-6)
})

test_that("a final minibatch of one row does not break training", {
  # 1801 training rows in batches of 200 leaves a single row over, and one
  # simulation has no contrasting parameter to be scored against. Before
  # minibatches() folded it in, backward() had no graph to walk and the fit
  # died with "element 0 of tensors does not require grad".
  skip_if_no_torch()
  sims <- simulate_for_sbi(gauss_sim, gauss_prior(), n = 2001L, seed = 1)
  expect_identical((2001L - 200L) %% 200L, 1L)

  fit <- nre(gauss_prior(), theta = sims$theta, x = sims$x, hidden = 8L,
             max_epochs = 2L, seed = 1)

  expect_true(is.finite(fit$de$best_val_loss))
})

test_that("the resnet classifier reproduces nflows' ResidualNet exactly", {
  # `sbi`'s "resnet" classifier is nflows' `ResidualNet`, which is what
  # nre_module() reimplements. Both were filled with the same deterministic
  # weights -- w[i, j] = sin(i + 2j + 1) / 3 for a matrix, b[i] = cos(i + 1) / 5
  # for a bias, zero-indexed -- and evaluated on the same three rows. The
  # reference below came out of
  # nflows.nn.nets.ResidualNet(in_features = 4, out_features = 1,
  #                            hidden_features = 6, num_blocks = 2,
  #                            activation = relu, dropout_probability = 0,
  #                            use_batch_norm = False)
  # under nflows 0.14. Anything that changes the block's shape -- where the
  # activations sit, which linear layer is zero-initialized, whether the final
  # layer sees an activation -- moves these numbers.
  skip_if_no_torch()
  net <- nre_module(2L, 2L, "resnet", hidden = 6L, n_blocks = 2L)()
  fill <- function(t) {
    d <- dim(t)
    torch::with_no_grad(t$copy_(torch::torch_tensor(
      if (length(d) == 2L) {
        sin(matrix(seq_len(d[1]) - 1L, d[1], d[2]) +
              2 * matrix(rep(seq_len(d[2]) - 1L, each = d[1]), d[1], d[2]) + 1) / 3
      } else {
        cos(seq_len(d[1])) / 5
      },
      dtype = torch::torch_float())))
  }
  for (p in net$parameters) fill(p)
  net$eval()

  z <- matrix(c(0.3, -0.4, 1.1, 0.0,
                -1.2, 0.7, 0.2, -0.9,
                2.0, 2.0, -2.0, 0.5), nrow = 3L, byrow = TRUE)
  got <- torch::with_no_grad(as.numeric(nre_logit_tensor(
    net,
    torch::torch_tensor(z[, 1:2], dtype = torch::torch_float()),
    torch::torch_tensor(z[, 3:4], dtype = torch::torch_float()))))

  expect_equal(got, c(0.108866699, 0.098798081, 0.180982038), tolerance = 1e-6)
})

test_that("a one-row batch scores zero through the graph, not around it", {
  # The guard the "final minibatch of one row" fix leans on. One simulation has
  # no contrast, so its atomic loss is zero -- but a zero torch did not build
  # carries no graph, and that is what made backward() fail. Returning it
  # through the classifier keeps the graph, and the gradient it produces is the
  # zero the batch is worth.
  skip_if_no_torch()
  net <- nre_module(1L, 1L, "mlp", hidden = 4L, n_blocks = 1L)()
  loss_fn <- nre_atomic_log_prob(10L)
  one <- torch::torch_tensor(matrix(0.5), dtype = torch::torch_float())

  lp <- loss_fn(net, one, one)

  expect_equal(as.numeric(lp), 0)
  expect_false(is.null(lp$grad_fn))
  expect_silent((-lp$mean())$backward())
  expect_true(all(as.numeric(net$parameters[[1]]$grad) == 0))
})
