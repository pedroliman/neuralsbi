test_that("named prior vectors attach param_names and propagate to sample_prior", {
  prior <- prior_uniform(low = c(beta = 0, gamma = 0), high = c(beta = 1, gamma = 1))
  expect_equal(prior$param_names, c("beta", "gamma"))
  th <- sample_prior(prior, 10)
  expect_equal(colnames(th), c("beta", "gamma"))

  prior2 <- prior_normal(mean = c(`beta[1]` = 0, rho = 0), sd = 1)
  expect_equal(prior2$param_names, c("beta[1]", "rho"))
  th2 <- sample_prior(prior2, 10)
  expect_equal(colnames(th2), c("beta[1]", "rho"))
})

test_that("as_theta_matrix preserves names for a single named vector row", {
  x <- c(alpha = 1, beta = 2)
  m <- as_theta_matrix(x, 2)
  expect_equal(colnames(m), c("alpha", "beta"))
  # a longer named vector reinterpreted as a stacked column loses its names
  x2 <- c(a = 1, b = 2, c = 3)
  m2 <- as_theta_matrix(x2, 1)
  expect_null(colnames(m2))
})

test_that("npe() carries theta/x names onto the fit, posterior samples, and MAP", {
  set.seed(1)
  prior <- prior_normal(mean = c(beta = 0, rho = 0), sd = 1)
  simulator <- function(theta) {
    out <- theta + matrix(rnorm(length(theta), sd = 0.3), nrow = nrow(theta))
    colnames(out) <- c("cases", "deaths")
    out
  }
  fit <- npe(prior, simulator, n_simulations = 1500,
             density_estimator = "linear_gaussian")
  expect_equal(fit$param_names, c("beta", "rho"))
  expect_equal(fit$x_names, c("cases", "deaths"))

  post <- posterior(fit, x_obs = c(cases = 0.5, deaths = -0.2))
  draws <- sample(post, 200)
  expect_equal(colnames(draws), c("beta", "rho"))

  map <- map_estimate(post, n_init = 200)
  expect_equal(names(map), c("beta", "rho"))

  pp <- posterior_predictive(post, simulator, n = 100)
  expect_equal(colnames(pp), c("cases", "deaths"))
})

test_that("sbc() and expected_coverage() carry parameter names", {
  set.seed(2)
  prior <- prior_normal(mean = c(beta = 0, rho = 0), sd = 1)
  simulator <- function(theta) theta + matrix(rnorm(length(theta), sd = 0.3),
                                              nrow = nrow(theta))
  fit <- npe(prior, simulator, n_simulations = 1500,
             density_estimator = "linear_gaussian")
  res <- sbc(fit, simulator, n_sbc = 50, n_posterior_samples = 100, seed = 3)
  expect_equal(colnames(res$ranks), c("beta", "rho"))
  expect_equal(names(res$uniformity_pvalue), c("beta", "rho"))

  cov <- expected_coverage(res, levels = c(0.5, 0.9))
  expect_equal(colnames(cov), c("nominal", "beta", "rho"))
})

test_that("math_expr parses valid R syntax and passes through otherwise", {
  expect_true(is.language(math_expr("beta[1]")))
  expect_true(is.language(math_expr("rho")))
  expect_identical(math_expr(NULL), NULL)
})

test_that("math_safe_text quotes unparseable labels and leaves valid ones alone", {
  out <- math_safe_text(c("beta[1]", "growth rate (per day)", "rho"))
  expect_equal(out[1], "beta[1]")
  expect_equal(out[3], "rho")
  # the unparseable label got quoted, and the quoted form must itself parse
  expect_true(grepl("^\".*\"$", out[2]))
  expect_silent(parse(text = out))
})

test_that("math_labels returns a same-length expression vector", {
  ex <- math_labels(c("beta[1]", "not valid !!", "rho"))
  expect_true(inherits(ex, "expression"))
  expect_length(ex, 3L)
})

test_that("plot_sbc, pairplot and plot_posterior_predictive run with named/math labels", {
  skip_if_no_ggally()
  set.seed(4)
  prior <- prior_normal(mean = c(`beta[1]` = 0, rho = 0), sd = 1)
  simulator <- function(theta) {
    out <- theta + matrix(rnorm(length(theta), sd = 0.3), nrow = nrow(theta))
    colnames(out) <- c("cases", "growth rate (per day)")
    out
  }
  fit <- npe(prior, simulator, n_simulations = 1500,
             density_estimator = "linear_gaussian")
  res <- sbc(fit, simulator, n_sbc = 30, n_posterior_samples = 100, seed = 5)

  path <- tempfile(fileext = ".png")
  grDevices::png(path)
  expect_silent(plot_sbc(res, param = 1L))
  post <- posterior(fit, x_obs = c(cases = 0.1, `growth rate (per day)` = -0.1))
  draws <- sample(post, 300)
  expect_silent(pairplot(draws))
  pp <- posterior_predictive(post, simulator, n = 200)
  expect_silent(plot_posterior_predictive(pp, c(0.1, -0.1)))
  grDevices::dev.off()
  expect_true(file.exists(path))
  unlink(path)
})

test_that("pairplot handles a non-syntactic label without erroring", {
  skip_if_no_ggally()
  set.seed(6)
  draws <- matrix(rnorm(200), ncol = 2)
  colnames(draws) <- c("growth rate (per day)", "rho")
  path <- tempfile(fileext = ".png")
  grDevices::png(path)
  expect_silent(pairplot(draws))
  grDevices::dev.off()
  expect_true(file.exists(path))
  unlink(path)
})
