# summary()/print()/as.data.frame() ergonomics and plot_coverage(), torch-free
# via the linear_gaussian estimator.

fit_lg <- function() {
  set.seed(7)
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + rnorm(length(theta), sd = 0.5)
  npe(prior, simulator, n_simulations = 1500,
      density_estimator = "linear_gaussian", seed = 7)
}

test_that("samples convert to data frames and summarize per parameter", {
  fit <- fit_lg()
  post <- posterior(fit, x_obs = c(1, -0.5))
  draws <- sample(post, 2000)

  df <- as.data.frame(draws)
  expect_s3_class(df, "data.frame")
  expect_equal(dim(df), c(2000L, 2L))
  expect_equal(names(df), c("theta1", "theta2"))

  s <- summary(draws)
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), 2L)
  expect_true(all(c("parameter", "mean", "sd", "q50") %in% names(s)))
  expect_true(all(s$q2.5 < s$q50 & s$q50 < s$q97.5))

  s2 <- summary(post, n = 500)
  expect_equal(nrow(s2), 2L)
})

test_that("summary of a fit returns training info invisibly", {
  fit <- fit_lg()
  info <- withVisible(summary(fit))
  expect_false(info$visible)
  expect_equal(info$value$density_estimator, "linear_gaussian")
  expect_equal(info$value$n_simulations, 1500L)
})

test_that("summary of an nle fit returns training info invisibly", {
  set.seed(11)
  prior <- prior_uniform(c(mu = -3, nu = -3), c(mu = 3, nu = 3))
  simulator <- function(mu, nu) c(a = mu + rnorm(1, sd = 0.4), b = nu + rnorm(1, sd = 0.3))
  fit <- nle(prior, simulator, n_simulations = 1200,
             density_estimator = "linear_gaussian", seed = 11)

  info <- withVisible(summary(fit))
  expect_false(info$visible)
  expect_equal(info$value$density_estimator, "linear_gaussian")
  expect_equal(info$value$n_simulations, 1200L)
  expect_equal(info$value$dim_theta, 2L)
  expect_equal(info$value$dim_x, 2L)
  expect_true(all(c("best_val_loss", "epochs_trained") %in% names(info$value)))
  expect_output(summary(fit), "Neural Likelihood Estimation")
})

test_that("printing a fit reports the estimator, the dimensions and the budget", {
  fit <- fit_lg()
  out <- paste(utils::capture.output(print(fit)), collapse = "\n")

  expect_match(out, "<nsbi_npe> Neural Posterior Estimation fit")
  expect_match(out, "density estimator : linear_gaussian")
  expect_match(out, "parameters \\(dim\\)  : 2")
  expect_match(out, "data \\(dim\\)        : 2")
  expect_match(out, "simulations       : 1500")
  expect_match(out, "posterior\\(fit, x_obs = \\.\\.\\.\\)")
  # The dead-network line belongs to a torch fit round-tripped through
  # saveRDS(), which test-serialize.R covers. A live fit must not print it.
  expect_false(grepl("network unusable", out))
})

test_that("printing a fit and its posterior lists the names when there are any", {
  set.seed(9)
  prior <- prior_uniform(c(mu = -3, nu = -3), c(mu = 3, nu = 3))
  simulator <- function(mu, nu) c(a = mu + rnorm(1, sd = 0.3), b = nu)
  fit <- npe(prior, simulator, n_simulations = 400,
             density_estimator = "linear_gaussian", seed = 9)

  name_lines <- grep("names", utils::capture.output(print(fit)), value = TRUE)
  expect_length(name_lines, 2L)
  expect_match(name_lines[1], "mu, nu")   # parameters
  expect_match(name_lines[2], "a, b")     # data

  out <- paste(utils::capture.output(print(posterior(fit, x_obs = c(0.5, -1)))),
               collapse = "\n")
  expect_match(out, "names         : mu, nu")
  expect_match(out, "conditioned on x: a=0.5, b=-1")
})

test_that("printing a posterior shows the observation it is conditioned on", {
  fit <- fit_lg()
  out <- paste(utils::capture.output(print(posterior(fit, x_obs = c(1, -0.5)))),
               collapse = "\n")

  expect_match(out, "<nsbi_posterior>")
  expect_match(out, "parameters \\(dim\\): 2")
  expect_match(out, "conditioned on x: 1, -0.5")
  expect_match(out, "sample\\(post, n\\), log_prob\\(post, theta\\)")

  # x_obs is optional, and the print has to say so rather than showing nothing.
  bare <- paste(utils::capture.output(print(posterior(fit))), collapse = "\n")
  expect_match(bare, "conditioned on x: \\(none set\\)")
})

test_that("printing samples reports the shape, the acceptance rate and a head", {
  fit <- fit_lg()
  draws <- sample(posterior(fit, x_obs = c(1, -0.5)), 500)
  lines <- utils::capture.output(print(draws))

  expect_match(lines[1], "<nsbi_samples> 500 draws x 2 parameters")
  # An unbounded prior rejects nothing, so every draw is kept.
  expect_match(lines[2], "support acceptance rate: 1\\.000")
  # The head prints as a plain matrix: the method must not recurse into itself.
  expect_gt(length(lines), 2L)
  expect_false(any(grepl("nsbi_samples", lines[-1])))
})

test_that("printing samples with no acceptance rate omits that line", {
  # An NLE posterior samples by MCMC rather than by rejection against the prior
  # support, so its draws carry no acceptance_rate attribute.
  set.seed(12)
  prior <- prior_uniform(c(mu = -3), c(mu = 3))
  fit <- nle(prior, function(mu) c(y = mu + rnorm(1, sd = 0.4)),
             n_simulations = 800, density_estimator = "linear_gaussian",
             seed = 12)
  post <- posterior(fit, matrix(0.5, nrow = 1), n_chains = 2, warmup = 50,
                    seed = 13)
  lines <- utils::capture.output(print(sample(post, 200)))

  expect_match(lines[1], "<nsbi_samples> 200 draws x 1 parameters")
  expect_false(any(grepl("acceptance rate", lines)))
})

test_that("plot_coverage runs on an sbc result and returns coverage", {
  skip_if_no_ggplot2()
  fit <- fit_lg()
  simulator <- function(theta) theta + rnorm(length(theta), sd = 0.5)
  res <- sbc(fit, simulator, n_sbc = 40L, n_posterior_samples = 200L, seed = 8)
  pdf(NULL)
  on.exit(dev.off())
  cov <- plot_coverage(res)
  expect_s3_class(cov, "data.frame")
  expect_equal(names(cov)[1], "nominal")
  # rough calibration for the exact estimator
  expect_lt(max(abs(cov$param1 - cov$nominal)), 0.3)
})
