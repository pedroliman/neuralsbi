test_that("uniform prior samples within bounds and has correct shape", {
  prior <- prior_uniform(low = c(-2, 0), high = c(2, 5))
  th <- sample_prior(prior, 100)
  expect_equal(dim(th), c(100L, 2L))
  expect_true(all(th[, 1] >= -2 & th[, 1] <= 2))
  expect_true(all(th[, 2] >= 0 & th[, 2] <= 5))
})

test_that("sample_prior() validates n", {
  prior <- prior_uniform(low = c(-2, 0), high = c(2, 5))
  expect_error(sample_prior(prior, 2.5), "`n` must be a single whole number")
  expect_error(sample_prior(prior, -1), "`n` must be a single whole number")
  expect_error(sample_prior(prior, NA), "`n` must be a single whole number")
})

test_that("uniform prior log_prob is constant inside, -Inf outside", {
  prior <- prior_uniform(low = c(0, 0), high = c(1, 1))
  inside <- matrix(c(0.5, 0.5), nrow = 1)
  outside <- matrix(c(1.5, 0.5), nrow = 1)
  expect_equal(prior$log_prob(inside), 0)          # log(1/(1*1)) = 0
  expect_equal(prior$log_prob(outside), -Inf)
})

test_that("normal prior log_prob matches dnorm", {
  prior <- prior_normal(mean = c(0, 1), sd = c(1, 2))
  th <- matrix(c(0, 1), nrow = 1)
  expect_equal(prior$log_prob(th),
               dnorm(0, 0, 1, log = TRUE) + dnorm(1, 1, 2, log = TRUE))
})

test_that("within_support flags out-of-bounds rows", {
  prior <- prior_uniform(c(-1, -1), c(1, 1))
  th <- rbind(c(0, 0), c(2, 0), c(0, -3))
  expect_equal(within_support(prior, th), c(TRUE, FALSE, FALSE))
})

test_that("unbounded prior treats everything as in-support", {
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  expect_true(all(within_support(prior, matrix(rnorm(20), ncol = 2))))
})

test_that("prior_uniform rejects bounds that do not describe an interval", {
  expect_error(prior_uniform(low = 2, high = 1),
               "Every `high` must be strictly greater than the matching `low`")
  expect_error(prior_uniform(low = c(0, 0), high = c(1, 0)),
               "Every `high` must be strictly greater")
  # An interval of zero width divides by zero in the log-density constant.
  expect_error(prior_uniform(low = 0, high = 0),
               "Every `high` must be strictly greater")
  expect_error(prior_uniform(low = c(0, 0), high = 1),
               "`low` and `high` must have the same length")
  # The same numbers the other way round are a valid prior.
  expect_s3_class(prior_uniform(low = 1, high = 2), "nsbi_prior")
})

test_that("prior_uniform rejects non-finite low/high but allows Inf", {
  expect_error(prior_uniform(low = NA, high = 5),
               "`low` contains 1 non-finite value \\(NA\\)")
  expect_error(prior_uniform(low = 0, high = NaN),
               "`high` contains 1 non-finite value \\(NaN\\)")
  expect_error(prior_uniform(low = c(0, NA), high = c(1, 2)),
               "`low` contains 1 non-finite value \\(NA\\)")
  expect_error(prior_uniform(low = NULL, high = 5), "`low` must be numeric")
  # An infinite bound is allowed at construction time, for prior_truncated()
  # to bound it.
  prior <- prior_uniform(low = c(mu = 0), high = c(mu = Inf))
  expect_s3_class(prior, "nsbi_prior")
  expect_equal(prior$upper, Inf)
})

test_that("sample_prior() errors on an improper (infinite-bound) uniform prior", {
  prior <- prior_uniform(low = c(mu = 0), high = c(mu = Inf))
  expect_error(sample_prior(prior, 5), "prior_truncated\\(\\)")
  # Fixed at both sides, sampling still works as always.
  finite <- prior_uniform(low = c(mu = 0), high = c(mu = 10))
  expect_s3_class(finite, "nsbi_prior")
  th <- sample_prior(finite, 5)
  expect_true(all(th >= 0 & th <= 10))
})

test_that("prior_normal rejects non-finite mean/sd, including Inf", {
  expect_error(prior_normal(mean = NA, sd = 1),
               "`mean` contains 1 non-finite value \\(NA\\)")
  expect_error(prior_normal(mean = 0, sd = NaN),
               "`sd` contains 1 non-finite value \\(NaN\\)")
  expect_error(prior_normal(mean = Inf, sd = 1),
               "`mean` contains 1 non-finite value \\(Inf\\)")
  expect_error(prior_normal(mean = 0, sd = Inf),
               "`sd` contains 1 non-finite value \\(Inf\\)")
  expect_error(prior_normal(mean = c(0, NA), sd = 1),
               "`mean` contains 1 non-finite value \\(NA\\)")
})

test_that("prior_normal rejects an sd that is not positive or not one per mean", {
  expect_error(prior_normal(mean = 0, sd = -1), "`sd` must be positive")
  expect_error(prior_normal(mean = 0, sd = 0), "`sd` must be positive")
  expect_error(prior_normal(mean = c(0, 0), sd = c(1, -2)),
               "`sd` must be positive")
  expect_error(prior_normal(mean = c(0, 0, 0), sd = c(1, 2)),
               "`sd` must be length 1 or the same length as `mean`")
  # A length-1 sd is recycled across every parameter, which is the case the
  # length check has to let through.
  expect_equal(prior_normal(mean = c(0, 0, 0), sd = 2)$params$sd, rep(2, 3))
})

test_that("prior_custom builds a working prior and keeps the names", {
  prior <- prior_custom(
    sample_fn = function(n) cbind(stats::runif(n), stats::runif(n, 0, 2)),
    log_prob_fn = function(theta) rep(-log(2), nrow(theta)),
    dim = 2L, lower = c(0, 0), upper = c(1, 2),
    param_names = c("beta", "gamma")
  )
  expect_equal(prior$param_names, c("beta", "gamma"))
  expect_equal(colnames(sample_prior(prior, 5)), c("beta", "gamma"))
  expect_equal(prior$log_prob(matrix(c(0.5, 1), nrow = 1)), -log(2))
})

test_that("prior_custom recycles a length-1 bound to one per parameter", {
  prior <- prior_custom(
    sample_fn = function(n) matrix(stats::runif(2 * n), ncol = 2),
    dim = 2L, lower = 0, upper = 1
  )
  expect_equal(prior$lower, c(0, 0))
  expect_equal(prior$upper, c(1, 1))
  expect_equal(within_support(prior, rbind(c(0.5, 0.5), c(0.5, 2))),
               c(TRUE, FALSE))
})

test_that("prior_custom rejects bounds that are not one per parameter", {
  sample_fn <- function(n) matrix(stats::runif(2 * n), ncol = 2)
  # sweep() used to recycle this and warn, leaving a wrong support test.
  expect_error(prior_custom(sample_fn, dim = 2L, lower = c(0, 0, 0)),
               "`lower` must be numeric of length 2")
  expect_error(prior_custom(sample_fn, dim = 2L, upper = "a"),
               "`upper` must be numeric of length 2")
  expect_error(prior_custom(sample_fn, dim = 2L, lower = c(0, NA)),
               "`lower` must be numeric of length 2")
  expect_error(
    prior_custom(sample_fn, dim = 2L, lower = c(0, 0), upper = c(1, 0)),
    "strictly greater"
  )
})

test_that("prior_custom rejects a sample_fn whose width disagrees with dim", {
  expect_error(
    prior_custom(function(n) matrix(stats::runif(n), ncol = 1), dim = 2L),
    "`sample_fn\\(2\\)` must return a 2 x 2 matrix"
  )
  # A sample_fn that ignores n returns the wrong number of rows.
  expect_error(
    prior_custom(function(n) matrix(stats::runif(6), ncol = 2), dim = 2L),
    "returned a 3 x 2 matrix"
  )
  expect_error(prior_custom(function(n) letters[seq_len(n)], dim = 1L),
               "`sample_fn\\(2\\)` must be numeric")
  expect_error(prior_custom(function(n) stop("boom"), dim = 1L),
               "`sample_fn` failed when called as `sample_fn\\(2\\)`: boom")
})

test_that("prior_custom rejects a log_prob_fn that is not one value per row", {
  sample_fn <- function(n) matrix(stats::runif(2 * n), ncol = 2)
  # A scalar return passed every probe downstream and surfaced as an MCMC
  # initialization failure, which blamed the sampler for a malformed prior.
  expect_error(
    prior_custom(sample_fn, log_prob_fn = function(theta) 0, dim = 2L),
    "must return one log-density per row of `theta`.*length-1 numeric vector"
  )
  expect_error(
    prior_custom(sample_fn, log_prob_fn = function(theta) "a", dim = 2L),
    "must return one log-density per row"
  )
  expect_error(
    prior_custom(sample_fn, log_prob_fn = function(theta) stop("nope"),
                 dim = 2L),
    "`log_prob_fn` failed on a 2 x 2 matrix"
  )
})

test_that("prior_custom checks dim, the callbacks and param_names", {
  sample_fn <- function(n) matrix(stats::runif(2 * n), ncol = 2)
  expect_error(prior_custom(sample_fn, dim = 0), "`dim` must be a single whole")
  expect_error(prior_custom(sample_fn, dim = 2.5), "`dim` must be a single whole")
  expect_error(prior_custom(sample_fn = "not a function", dim = 2L),
               "`sample_fn` must be a function of one argument")
  expect_error(prior_custom(function() 1, dim = 2L),
               "`sample_fn` must be a function of one argument .*takes none")
  expect_error(prior_custom(sample_fn, log_prob_fn = 1, dim = 2L),
               "`log_prob_fn` must be a function of one argument")
  expect_error(prior_custom(sample_fn, dim = 2L, param_names = "beta"),
               "`param_names` must be a character vector with one non-empty")
  expect_error(prior_custom(sample_fn, dim = 2L, param_names = c("beta", "")),
               "`param_names` must be a character vector")
})

test_that("printing a one-sided custom prior shows the bound it has", {
  prior <- prior_custom(function(n) matrix(stats::rexp(2 * n), ncol = 2),
                        dim = 2L, lower = 0)
  out <- paste(utils::capture.output(print(prior)), collapse = "\n")
  expect_match(out, "lower: 0, 0")
  expect_false(grepl("upper", out))
})

test_that("a custom prior without log_prob_fn still returns NA per row", {
  prior <- prior_custom(function(n) matrix(stats::runif(2 * n), ncol = 2),
                        dim = 2L)
  expect_equal(prior$log_prob(matrix(0, nrow = 3, ncol = 2)),
               rep(NA_real_, 3))
})
