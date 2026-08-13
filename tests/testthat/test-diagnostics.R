test_that("sbc() rejects an object that is not an npe/nle/nre fit, naming its class", {
  expect_error(
    sbc(structure(list(), class = "lm"), function(theta) theta,
        prior_normal(mean = 0, sd = 1)),
    "Expected a fit from npe\\(\\), nle\\(\\) or nre\\(\\), not an object of class lm"
  )
})

test_that("sbc returns ranks of the right shape and reasonable calibration", {
  set.seed(7)
  d <- 2; sigma <- 0.5
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + rnorm(length(theta), sd = sigma)
  fit <- npe(prior, simulator, n_simulations = 4000,
             density_estimator = "linear_gaussian")
  res <- sbc(fit, simulator, n_sbc = 100, n_posterior_samples = 200, seed = 3)
  expect_equal(dim(res$ranks), c(100L, 2L))
  # a well-specified exact estimator should not fail uniformity badly
  expect_true(all(res$uniformity_pvalue > 0.001))
})

test_that("print.nsbi_sbc() labels each parameter's p-value when the fit is named", {
  set.seed(2)
  prior <- prior_normal(mean = c(beta = 0, rho = 0), sd = 1)
  simulator <- function(theta) unname(theta) + rnorm(2, sd = 0.3)
  fit <- npe(prior, simulator, n_simulations = 1500,
             density_estimator = "linear_gaussian")
  res <- sbc(fit, simulator, n_sbc = 50, n_posterior_samples = 100, seed = 3)
  expect_output(print(res), "beta=")
  expect_output(print(res), "rho=")
})

test_that("expected_coverage produces a monotone-ish curve near the diagonal", {
  set.seed(7)
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + rnorm(length(theta), sd = 0.5)
  fit <- npe(prior, simulator, n_simulations = 4000,
             density_estimator = "linear_gaussian")
  res <- sbc(fit, simulator, n_sbc = 200, n_posterior_samples = 200, seed = 5)
  cov <- expected_coverage(res, levels = c(0.5, 0.9))
  expect_true(all(cov$param1 >= 0 & cov$param1 <= 1))
  # 90% interval should cover clearly more often than the 50% interval
  expect_gt(cov$param1[2], cov$param1[1])
})

test_that("expected_coverage() refuses levels outside (0, 1)", {
  # A level outside the unit interval used to be scored anyway: the central
  # interval comes out empty or covers the line, so the row reads as coverage
  # 0 or 1 and looks like a verdict on the fit.
  set.seed(14)
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + stats::rnorm(length(theta), sd = 0.5)
  fit <- npe(prior, simulator, n_simulations = 400,
             density_estimator = "linear_gaussian")
  res <- sbc(fit, simulator, n_sbc = 100L, n_posterior_samples = 50L, seed = 2)

  expect_error(expected_coverage(res, levels = c(-1, 2)),
               regexp = "`levels` must be numbers strictly between 0 and 1")
  # The values are listed, so which entry is wrong is readable.
  expect_error(expected_coverage(res, levels = c(0.5, 2)), regexp = "0.5, 2")
  expect_error(expected_coverage(res, levels = c(0.5, NA)), regexp = "0.5, NA")
  expect_error(expected_coverage(res, levels = 0), regexp = "not 0")
  expect_error(expected_coverage(res, levels = 1), regexp = "not 1")
  expect_error(expected_coverage(res, levels = "0.5"),
               regexp = "not a character value")
  expect_error(expected_coverage(res, levels = numeric(0)),
               regexp = "not a length-0 numeric vector")

  cov <- expected_coverage(res, levels = c(0.5, 0.9))
  expect_equal(cov$nominal, c(0.5, 0.9))
  expect_true(all(cov$param1 >= 0 & cov$param1 <= 1))
})

test_that("sbc errors instead of ranking against draws it never got", {
  # A bounded prior plus an estimator that leaks means sample() comes back
  # short, and sbc() used to bin those ranks against n_posterior_samples: every
  # rank compressed toward zero and the run read as miscalibrated. Force the
  # leak by pushing the fitted conditional mean 50 standardized units off the
  # prior box, so every draw is rejected.
  set.seed(11)
  prior <- prior_uniform(low = c(0, 0), high = c(1, 1))
  simulator <- function(theta) theta + rnorm(length(theta), sd = 0.2)
  fit <- npe(prior, simulator, n_simulations = 500,
             density_estimator = "linear_gaussian")
  fit$de$B[1, ] <- fit$de$B[1, ] + 50

  expect_error(
    suppressWarnings(sbc(fit, simulator, n_sbc = 2L,
                         n_posterior_samples = 20L, seed = 1)),
    "0 of 20 posterior draws"
  )
})

test_that("tarp errors on a short posterior draw too", {
  set.seed(12)
  prior <- prior_uniform(low = c(0, 0), high = c(1, 1))
  simulator <- function(theta) theta + rnorm(length(theta), sd = 0.2)
  fit <- npe(prior, simulator, n_simulations = 500,
             density_estimator = "linear_gaussian")
  fit$de$B[1, ] <- fit$de$B[1, ] + 50

  expect_error(
    suppressWarnings(tarp(fit, simulator, n_tarp = 2L,
                          n_posterior_samples = 20L, seed = 1)),
    "posterior draws"
  )
})

test_that("c2st of a sample set against itself is ~0.5", {
  set.seed(1)
  a <- matrix(rnorm(2000), ncol = 2)
  b <- matrix(rnorm(2000), ncol = 2)
  res <- c2st(a, b, seed = 1)
  expect_lt(res$accuracy, 0.6)
})

test_that("c2st is not fooled by unequal sample sizes", {
  # Accuracy against unbalanced classes is not a two-sample test. Four times as
  # many draws on one side and a classifier scores 0.8 by always answering with
  # the bigger class, having learned nothing -- which is what the vignette hit
  # comparing 8000 slice draws against 2000 from Stan.
  set.seed(2)
  a <- matrix(rnorm(16000), ncol = 2)
  b <- matrix(rnorm(4000), ncol = 2)

  expect_lt(c2st(a, b, seed = 1)$accuracy, 0.6)
  expect_lt(c2st(b, a, seed = 1)$accuracy, 0.6)

  # Balancing must not cost it the ability to see a real difference.
  shifted <- matrix(rnorm(4000, mean = 3), ncol = 2)
  expect_gt(c2st(a, shifted, seed = 1)$accuracy, 0.9)
})

test_that("posterior_predictive returns simulator-shaped output", {
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + rnorm(1, sd = 0.3)
  fit <- npe(prior, simulator, n_simulations = 1000,
             density_estimator = "linear_gaussian")
  post <- posterior(fit, x_obs = 0.5)
  pp <- posterior_predictive(post, simulator, n = 300)
  expect_equal(nrow(pp), 300L)
})

test_that("posterior_predictive() errors clearly when the posterior has no draws in support", {
  # An observation far outside anything the bounded prior can explain drives
  # rejection sampling to zero draws (issue #171). Before the fix this hit
  # `dimnames<-` on a 0x0 matrix and errored with a message that named neither
  # the function nor the actual problem.
  set.seed(1)
  prior <- prior_uniform(c(mu = -1, nu = -1), c(mu = 1, nu = 1))
  simulator <- function(mu, nu) c(a = mu + rnorm(1, sd = 0.05), b = nu + rnorm(1, sd = 0.05))
  fit <- npe(prior, simulator, n_simulations = 500,
             density_estimator = "linear_gaussian", seed = 1)
  post <- posterior(fit, x_obs = c(6, 6))
  expect_error(
    suppressWarnings(posterior_predictive(post, simulator, n = 50, x = c(6, 6))),
    "posterior_predictive\\(\\): the posterior returned no draws"
  )
})

test_that("plot_posterior_predictive runs and locates the observation", {
  skip_if_no_ggplot2()
  set.seed(2)
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + rnorm(length(theta), sd = 0.3)
  fit <- npe(prior, simulator, n_simulations = 1000,
             density_estimator = "linear_gaussian")
  x_obs <- c(0.5, -0.2)
  post <- posterior(fit, x_obs = x_obs)
  pp <- posterior_predictive(post, simulator, n = 500)
  path <- tempfile(fileext = ".png")
  grDevices::png(path)
  q <- plot_posterior_predictive(pp, x_obs)
  grDevices::dev.off()
  expect_true(file.exists(path))
  expect_length(q, 2L)
  # the observation should sit inside the bulk of its own predictive
  expect_true(all(q > 0.01 & q < 0.99))
  unlink(path)
})

test_that("sbc(), tarp() and c2st() check their counts", {
  set.seed(11)
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta + stats::rnorm(1, sd = 0.5)
  fit <- npe(prior, simulator, n_simulations = 300,
             density_estimator = "linear_gaussian")

  expect_error(sbc(fit, simulator, n_sbc = 0), "`n_sbc` must be")
  expect_error(sbc(fit, simulator, n_sbc = 10, n_posterior_samples = 0),
               "`n_posterior_samples` must be")
  expect_error(tarp(fit, simulator, n_tarp = -1), "`n_tarp` must be")
  expect_error(tarp(fit, simulator, n_tarp = 10, n_posterior_samples = 2.5),
               "`n_posterior_samples` must be")

  draws <- matrix(stats::rnorm(100), ncol = 2)
  expect_error(c2st(draws, draws, n_folds = 1),
               "`n_folds` must be a single whole number of at least 2 since")
  expect_error(c2st(draws, draws, n_folds = NA), "`n_folds` must be")
})

test_that("c2st() refuses two sample sets of different widths", {
  set.seed(3)
  a <- matrix(stats::rnorm(200), ncol = 2)
  b <- matrix(stats::rnorm(300), ncol = 3)

  expect_error(c2st(a, b), regexp = "`x` has 2 columns and `y` has 3 columns")
  expect_error(c2st(b, a), regexp = "`x` has 3 columns and `y` has 2 columns")
  expect_error(c2st(a, letters[1:4]), regexp = "`y` must be numeric")
})

test_that("c2st() refuses more folds than it has draws to fill them", {
  # Left alone, rep_len() leaves the last folds empty, mean() of an empty test
  # fold is NaN, and the accuracy comes back NaN with nothing said about why.
  set.seed(4)
  few <- matrix(stats::rnorm(8), ncol = 2)

  expect_error(c2st(few, few), regexp = "smaller sample set has only 4 draws")
  expect_error(c2st(few, few, n_folds = 4),
               regexp = "`n_folds` is 4, but the smaller sample set")
  # The lower bound still comes from check_count(), before the draws are seen.
  expect_error(c2st(few, few, n_folds = 1),
               regexp = "at least 2 since each fold is scored")

  # Fewer folds than draws is fine, including on the smaller of two sets.
  expect_type(c2st(few, few, n_folds = 3, seed = 1)$accuracy, "double")
  plenty <- matrix(stats::rnorm(400), ncol = 2)
  expect_false(is.nan(c2st(plenty, few, n_folds = 3, seed = 1)$accuracy))
})

test_that("sbc() and tarp() refuse a prior that is not the width of the fit", {
  # `prior` exists to be overridden, so the width it comes with has to be
  # checked. A narrower one used to reach sweep() in sbc() and the z-scoring in
  # tarp(), where it recycles against the fit's width and the diagnostic scores
  # a comparison nobody asked for.
  set.seed(13)
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + stats::rnorm(length(theta), sd = 0.5)
  fit <- npe(prior, simulator, n_simulations = 400,
             density_estimator = "linear_gaussian")

  wrong <- prior_normal(mean = c(0, 0, 0), sd = 1)
  expect_error(sbc(fit, simulator, prior = wrong, n_sbc = 2L),
               regexp = "`prior` covers 3 parameters but 2 are expected here")
  expect_error(tarp(fit, simulator, prior = wrong, n_tarp = 2L),
               regexp = "`prior` covers 3 parameters but 2 are expected here")

  expect_error(sbc(fit, simulator, prior = prior_normal(mean = 0, sd = 1),
                   n_sbc = 2L),
               regexp = "covers 1 parameter but 2 are expected")
  expect_error(sbc(fit, simulator, prior = list(dim = 2L), n_sbc = 2L),
               regexp = "`prior` must be an nsbi_prior object, not list")

  # The check runs before the simulator does, since that is the expensive part.
  calls <- 0L
  counting <- function(theta) {
    calls <<- calls + 1L
    theta + stats::rnorm(length(theta), sd = 0.5)
  }
  expect_error(sbc(fit, counting, prior = wrong, n_sbc = 50L))
  expect_error(tarp(fit, counting, prior = wrong, n_tarp = 50L))
  expect_identical(calls, 0L)

  # A prior supplied explicitly at the right width still works, so the check
  # does not simply reject every override.
  narrower <- prior_normal(mean = c(0, 0), sd = 0.5)
  # 100 trials over the 20 rank bins keeps chisq.test() off its
  # small-expected-count warning; the point here is the shape, not the p-value.
  res <- sbc(fit, simulator, prior = narrower, n_sbc = 100L,
             n_posterior_samples = 50L, seed = 2)
  expect_equal(dim(res$ranks), c(100L, 2L))
  tp <- tarp(fit, simulator, prior = narrower, n_tarp = 20L,
             n_posterior_samples = 50L, seed = 2)
  expect_s3_class(tp, "nsbi_tarp")
})
