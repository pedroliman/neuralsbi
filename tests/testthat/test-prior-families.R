test_that("every family samples inside its support with the right shape", {
  priors <- list(
    lognormal = prior_lognormal(c(0, 1), c(0.5, 0.2)),
    exponential = prior_exponential(c(1, 3)),
    gamma = prior_gamma(c(2, 9), c(1, 3)),
    beta = prior_beta(c(2, 5), c(15, 1)),
    student_t = prior_student_t(3, c(0, 1), c(1, 2)),
    cauchy = prior_cauchy(c(0, 1), c(1, 2)),
    half_normal = prior_half_normal(c(1, 2)),
    half_cauchy = prior_half_cauchy(c(1, 2))
  )
  for (nm in names(priors)) {
    prior <- priors[[nm]]
    expect_identical(prior$type, nm)
    expect_identical(prior$dim, 2L)
    theta <- sample_prior(prior, 200)
    expect_equal(dim(theta), c(200L, 2L), info = nm)
    expect_true(all(is.finite(theta)), info = nm)
    expect_true(all(within_support(prior, theta)), info = nm)
  }
})

test_that("family log_prob matches the base R density", {
  theta <- cbind(c(0.4, 1.2, 3.0), c(0.1, 0.7, 2.5))

  expect_equal(
    prior_lognormal(c(0, 1), c(0.5, 0.2))$log_prob(theta),
    stats::dlnorm(theta[, 1], 0, 0.5, log = TRUE) +
      stats::dlnorm(theta[, 2], 1, 0.2, log = TRUE)
  )
  expect_equal(
    prior_exponential(c(1, 3))$log_prob(theta),
    stats::dexp(theta[, 1], 1, log = TRUE) +
      stats::dexp(theta[, 2], 3, log = TRUE)
  )
  expect_equal(
    prior_gamma(c(2, 9), c(1, 3))$log_prob(theta),
    stats::dgamma(theta[, 1], shape = 2, rate = 1, log = TRUE) +
      stats::dgamma(theta[, 2], shape = 9, rate = 3, log = TRUE)
  )
  expect_equal(
    prior_cauchy(c(0, 1), c(1, 2))$log_prob(theta),
    stats::dcauchy(theta[, 1], 0, 1, log = TRUE) +
      stats::dcauchy(theta[, 2], 1, 2, log = TRUE)
  )

  # Stan's student_t is the location-scale one; stats::dt() is not.
  expect_equal(
    prior_student_t(3, c(0, 1), c(1, 2))$log_prob(theta),
    stats::dt(theta[, 1], df = 3, log = TRUE) +
      stats::dt((theta[, 2] - 1) / 2, df = 3, log = TRUE) - log(2)
  )

  unit <- cbind(c(0.1, 0.5, 0.9), c(0.2, 0.4, 0.6))
  expect_equal(
    prior_beta(2, 15)$log_prob(unit[, 1, drop = FALSE]),
    stats::dbeta(unit[, 1], 2, 15, log = TRUE)
  )
})

test_that("the half families carry the log(2) renormalization", {
  x <- c(0.25, 1.5, 3)
  expect_equal(prior_half_normal(2)$log_prob(x),
               stats::dnorm(x, 0, 2, log = TRUE) + log(2))
  expect_equal(prior_half_cauchy(1.5)$log_prob(x),
               stats::dcauchy(x, 0, 1.5, log = TRUE) + log(2))

  # And they are proper densities on [0, Inf).
  for (prior in list(prior_half_normal(2), prior_half_cauchy(1.5))) {
    mass <- stats::integrate(function(z) exp(prior$log_prob(z)), 0, Inf)
    expect_equal(mass$value, 1, tolerance = 1e-5)
  }
})

test_that("a scalar argument is shared and the longest one sets the dimension", {
  prior <- prior_gamma(shape = c(2, 5, 9), rate = 3)
  expect_identical(prior$dim, 3L)
  rates <- vapply(prior$params$marginals, function(m) m$args$rate, numeric(1))
  expect_equal(rates, rep(3, 3))

  expect_error(prior_gamma(shape = c(2, 5), rate = c(1, 2, 3)),
               "Every argument to prior_gamma\\(\\) must be length 1 or the ")
  expect_error(prior_student_t(df = c(3, 4), location = c(0, 1, 2),
                               scale = c(1, 2, 3, 4)),
               "\\(4\\), but `df` and `location` are not")
})

test_that("bounded families set the support that drives leakage correction", {
  expect_equal(prior_lognormal(0, 1)$lower, 0)
  expect_null(prior_lognormal(0, 1)$upper)
  expect_equal(prior_beta(2, 15)$lower, 0)
  expect_equal(prior_beta(2, 15)$upper, 1)
  expect_equal(prior_half_normal(1)$lower, 0)
  expect_null(prior_cauchy(0, 1)$lower)
  expect_null(prior_cauchy(0, 1)$upper)

  expect_equal(within_support(prior_beta(2, 15), c(0.1, 1.5, -0.2)),
               c(TRUE, FALSE, FALSE))
  expect_equal(prior_beta(2, 15)$log_prob(c(0.1, 1.5))[2], -Inf)
})

test_that("parameter names come off whichever argument carries them", {
  prior <- prior_lognormal(meanlog = c(beta = log(0.4), gamma = log(0.125)),
                           sdlog = c(0.5, 0.2))
  expect_equal(prior$param_names, c("beta", "gamma"))
  expect_equal(colnames(sample_prior(prior, 4)), c("beta", "gamma"))

  # The first complete set wins, and a recycled scalar is not a complete set.
  expect_equal(prior_gamma(shape = c(2, 3), rate = c(a = 1))$param_names, NULL)
  expect_equal(prior_gamma(shape = c(2, 3), rate = c(a = 1, b = 2))$param_names,
               c("a", "b"))
})

test_that("family constructors reject parameters the d/q functions would not", {
  expect_error(prior_gamma(shape = 0), "`shape` must be one or more finite ")
  expect_error(prior_gamma(shape = c(2, -1)),
               "`shape` must be one or more finite positive numbers, not 2, -1")
  expect_error(prior_beta(2, NA), "`shape2` must be one or more finite ")
  expect_error(prior_lognormal(meanlog = Inf),
               "`meanlog` must be one or more finite numbers, not Inf")
  expect_error(prior_lognormal(sdlog = 0), "`sdlog` must be")
  expect_error(prior_exponential("fast"), "`rate` must be")
  expect_error(prior_student_t(df = 3, scale = -1), "`scale` must be")
  expect_error(prior_half_normal(sd = 0), "`sd` must be")
  expect_error(prior_half_cauchy(scale = numeric(0)), "`scale` must be")
})

# ---- prior_independent() ---------------------------------------------------

test_that("prior_independent multiplies the marginals and stacks the bounds", {
  prior <- prior_independent(
    p = prior_beta(2, 15),
    hr = prior_lognormal(log(3), 0.3),
    z = prior_normal(mean = 0, sd = 2)
  )
  expect_identical(prior$type, "independent")
  expect_identical(prior$dim, 3L)
  expect_equal(prior$param_names, c("p", "hr", "z"))
  expect_equal(prior$lower, c(0, 0, -Inf))
  expect_equal(prior$upper, c(1, Inf, Inf))

  theta <- cbind(c(0.1, 0.3), c(2, 5), c(-1, 1))
  expect_equal(
    prior$log_prob(theta),
    stats::dbeta(theta[, 1], 2, 15, log = TRUE) +
      stats::dlnorm(theta[, 2], log(3), 0.3, log = TRUE) +
      stats::dnorm(theta[, 3], 0, 2, log = TRUE)
  )
  draws <- sample_prior(prior, 100)
  expect_equal(colnames(draws), c("p", "hr", "z"))
  expect_true(all(within_support(prior, draws)))
})

test_that("prior_independent names a multi-parameter component", {
  # The component's own names win where it has them.
  named <- prior_independent(prior_normal(mean = c(a = 0, b = 0)),
                             s = prior_half_normal(1))
  expect_equal(named$param_names, c("a", "b", "s"))

  # Otherwise the argument name is numbered.
  numbered <- prior_independent(mu = prior_normal(mean = c(0, 0)),
                                s = prior_half_normal(1))
  expect_equal(numbered$param_names, c("mu1", "mu2", "s"))

  # A partial set is no set: a prior cannot label some columns and not others.
  expect_null(prior_independent(prior_normal(mean = c(0, 0)),
                                s = prior_half_normal(1))$param_names)
})

test_that("prior_independent composes a prior_custom component too", {
  custom <- prior_custom(
    sample_fn = function(n) matrix(stats::runif(n, 2, 4), ncol = 1),
    log_prob_fn = function(theta) rep(-log(2), nrow(theta)),
    dim = 1L, lower = 2, upper = 4
  )
  prior <- prior_independent(a = prior_beta(2, 15), b = custom)

  expect_identical(prior$dim, 2L)
  expect_equal(prior$lower, c(0, 2))
  expect_equal(prior$upper, c(1, 4))
  theta <- cbind(c(0.1, 0.4), c(3, 3.5))
  expect_equal(prior$log_prob(theta),
               stats::dbeta(theta[, 1], 2, 15, log = TRUE) - log(2))
  draws <- sample_prior(prior, 50)
  expect_equal(colnames(draws), c("a", "b"))
  expect_true(all(within_support(prior, draws)))

  # No family behind the custom half, so nothing to renormalize against.
  expect_null(prior$params$marginals)
  expect_error(prior_truncated(prior, upper = 0.5), "no family to renormalize")
})

test_that("prior_independent checks that it was given priors", {
  expect_error(prior_independent(), "needs at least one prior")
  expect_error(prior_independent(prior_beta(2, 2), 3),
               "argument 2 is of class numeric")
  expect_error(prior_independent(a = prior_beta(2, 2), b = NULL),
               "`b` is of class NULL")
})

# ---- prior_truncated() -----------------------------------------------------

test_that("prior_truncated renormalizes to a proper density", {
  prior <- prior_truncated(prior_lognormal(log(0.4), 0.5),
                           lower = 0.1, upper = 2)
  expect_identical(prior$type, "truncated")
  expect_equal(prior$lower, 0.1)
  expect_equal(prior$upper, 2)

  mass <- stats::integrate(function(z) exp(prior$log_prob(z)), 0.1, 2)
  expect_equal(mass$value, 1, tolerance = 1e-6)

  z <- log(0.4)
  const <- log(stats::plnorm(2, z, 0.5) - stats::plnorm(0.1, z, 0.5))
  expect_equal(prior$log_prob(c(0.5, 1)),
               stats::dlnorm(c(0.5, 1), z, 0.5, log = TRUE) - const)
  expect_equal(prior$log_prob(c(0.05, 3)), c(-Inf, -Inf))

  draws <- sample_prior(prior, 500)
  expect_true(all(draws >= 0.1 & draws <= 2))
})

test_that("a truncated normal is the half-normal, and truncation composes", {
  half <- prior_truncated(prior_normal(mean = 0, sd = 1), lower = 0)
  expect_equal(half$log_prob(c(0.5, 2)),
               prior_half_normal(1)$log_prob(c(0.5, 2)))

  # Truncating twice intersects the bounds rather than stacking constants.
  twice <- prior_truncated(prior_truncated(prior_normal(0, 1), lower = 0),
                           upper = 1)
  once <- prior_truncated(prior_normal(0, 1), lower = 0, upper = 1)
  expect_equal(twice$log_prob(0.5), once$log_prob(0.5))
  mass <- stats::integrate(function(z) exp(twice$log_prob(z)), 0, 1)
  expect_equal(mass$value, 1, tolerance = 1e-6)
})

test_that("prior_truncated bounds each parameter separately", {
  prior <- prior_truncated(prior_normal(mean = c(a = 0, b = 0), sd = 1),
                           lower = c(-1, 0), upper = c(1, Inf))
  expect_equal(prior$param_names, c("a", "b"))
  expect_equal(prior$lower, c(-1, 0))
  expect_equal(prior$upper, c(1, Inf))
  draws <- sample_prior(prior, 300)
  expect_true(all(draws[, 1] >= -1 & draws[, 1] <= 1))
  expect_true(all(draws[, 2] >= 0))
})

test_that("prior_truncated checks its arguments and its bounds", {
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  expect_error(prior_truncated(prior), "Give `lower`, `upper` or both")
  expect_error(prior_truncated(prior, lower = c(0, 0, 0)),
               "`lower` must be numeric of length 2")
  expect_error(prior_truncated(prior, lower = 1, upper = 0),
               "Every `upper` must be strictly greater")
  expect_error(prior_truncated(list(), lower = 0),
               "must be an nsbi_prior object")

  # Bounds off in the far tail leave no mass to renormalize by.
  expect_error(prior_truncated(prior_normal(0, 1), lower = 100, upper = 200),
               "leave no probability mass")
  # And bounds outside the family's own support leave nothing at all.
  expect_error(prior_truncated(prior_beta(2, 2), lower = 2, upper = 3),
               "is empty once it is intersected")

  custom <- prior_custom(function(n) matrix(stats::runif(n), ncol = 1),
                         dim = 1L, lower = 0, upper = 1)
  expect_error(prior_truncated(custom, upper = 0.5),
               "prior_truncated\\(\\) needs a prior built from a named family")
  expect_error(prior_truncated(custom, upper = 0.5), "prior_custom\\(\\)")
})

test_that("prior_uniform and prior_normal can be truncated", {
  half_box <- prior_truncated(prior_uniform(0, 4), upper = 1)
  expect_equal(half_box$log_prob(c(0.5, 2)), c(0, -Inf))
  expect_true(all(sample_prior(half_box, 100) <= 1))

  # An improper uniform has no CDF to invert, so it carries no family.
  expect_error(prior_truncated(prior_uniform(0, Inf), upper = 1),
               "no family to renormalize")
})

# ---- printing --------------------------------------------------------------

test_that("printing a family prior shows one line per marginal", {
  prior <- prior_independent(p = prior_beta(2, 15),
                             s = prior_half_normal(1.5))
  out <- paste(utils::capture.output(print(prior)), collapse = "\n")
  expect_match(out, "type=independent, dim=2")
  expect_match(out, "p ~ beta\\(2, 15\\)")
  expect_match(out, "s ~ normal\\(0, 1.5\\) T\\[0, \\]")

  unnamed <- utils::capture.output(print(prior_exponential(c(1, 2))))
  expect_match(paste(unnamed, collapse = "\n"), "theta\\[2\\] ~ exponential\\(2\\)")
})

test_that("prior_uniform and prior_normal print as they always have", {
  box <- utils::capture.output(print(prior_uniform(c(a = 0), c(a = 1))))
  expect_equal(box, c("<nsbi_prior> type=uniform, dim=1",
                      "  parameters: a ", "  lower: 0 ", "  upper: 1 "))
  norm <- utils::capture.output(print(prior_normal(mean = 0, sd = 1)))
  expect_equal(norm, "<nsbi_prior> type=normal, dim=1")
})

# ---- Stan emission ---------------------------------------------------------

test_that("stan_prior_blocks() writes a family out as a sampling statement", {
  blocks <- stan_prior_blocks(
    prior_independent(prior_lognormal(log(0.4), 0.5), prior_half_normal(2)), 2L)

  expect_equal(blocks$data, "")
  expect_equal(blocks$transformed, "")
  expect_match(blocks$parameters, "vector<lower=0>\\[2\\] theta;")
  expect_match(blocks$model, "theta\\[1\\] ~ lognormal\\(-0.916290731874155, 0.5\\);")
  expect_match(blocks$model, "theta\\[2\\] ~ normal\\(0, 2\\) T\\[0,\\];")
})

test_that("bounds that differ per parameter are declared one at a time", {
  blocks <- stan_prior_blocks(
    prior_independent(prior_gamma(2, 3), prior_cauchy(0, 2)), 2L)

  expect_match(blocks$parameters, "real<lower=0> theta_1;")
  expect_match(blocks$parameters, "real theta_2;")
  expect_match(blocks$transformed,
               "vector\\[2\\] theta = \\[theta_1, theta_2\\]';")
  expect_match(blocks$model, "theta_1 ~ gamma\\(2, 3\\);")
  expect_match(blocks$model, "theta_2 ~ cauchy\\(0, 2\\);")
})

test_that("truncation is written out as T[,] on the side that was cut", {
  blocks <- stan_prior_blocks(
    prior_truncated(prior_beta(2, 15), upper = 0.5), 1L)
  expect_match(blocks$parameters, "vector<lower=0, upper=0.5>\\[1\\] theta;")
  expect_match(blocks$model, "theta\\[1\\] ~ beta\\(2, 15\\) T\\[,0.5\\];")

  both <- stan_prior_blocks(
    prior_truncated(prior_student_t(3, 0, 1), lower = -2, upper = 2), 1L)
  expect_match(both$model,
               "theta\\[1\\] ~ student_t\\(3, 0, 1\\) T\\[-2,2\\];")
})

test_that("a truncated uniform stays implicit in its declared bounds", {
  blocks <- stan_prior_blocks(prior_truncated(prior_uniform(0, 4), upper = 1),
                              1L)
  expect_match(blocks$parameters, "vector<lower=0, upper=1>\\[1\\] theta;")
  expect_equal(blocks$model, "")
})

test_that("prior_uniform and prior_normal still travel through the data block", {
  box <- stan_prior_blocks(prior_uniform(c(0, 0), c(1, 1)), 2L)
  expect_match(box$data, "vector\\[2\\] nsbi_low;")
  expect_match(box$parameters, "vector<lower=nsbi_low, upper=nsbi_high>")
  expect_equal(box$model, "")

  norm <- stan_prior_blocks(prior_normal(c(0, 0), 1), 2L)
  expect_match(norm$data, "nsbi_prior_mean")
  expect_equal(norm$model, "  theta ~ normal(nsbi_prior_mean, nsbi_prior_sd);\n")
})

test_that("a custom prior still says to write the model block yourself", {
  custom <- prior_custom(function(n) matrix(stats::runif(2 * n), ncol = 2),
                         dim = 2L)
  expect_error(stan_prior_blocks(custom, 2L),
               "write the model block yourself")
  expect_error(stan_prior_blocks(custom, 2L), "\\?prior_families")
})

test_that("stan_code() puts the generated prior into a whole program", {
  fit <- nle(prior_independent(mu = prior_normal(mean = 0, sd = 2),
                               s = prior_half_normal(1)),
             function(theta) c(y = stats::rnorm(1, theta[1], 0.5),
                               z = stats::rnorm(1, theta[2], 0.5)),
             n_simulations = 300, density_estimator = "linear_gaussian",
             seed = 1)
  code <- stan_code(fit)

  expect_match(code, "transformed parameters \\{")
  expect_match(code, "theta_1 ~ normal\\(0, 2\\);")
  expect_match(code, "theta_2 ~ normal\\(0, 1\\) T\\[0,\\];")
  expect_match(code, "x ~ nsbi_log_lik_sum\\(theta, nsbi_w\\);")
  # Nothing extra to ship: the prior is in the source, not in the data.
  expect_named(stan_data(fit), c("nsbi_nw", "nsbi_w"))
})
