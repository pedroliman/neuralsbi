# save_npe()/load_npe(): the round trip that survives a torch back end.

fit_toy <- function(n = 500) {
  prior <- prior_normal(mean = c(mu = 0, nu = 0), sd = 1)
  npe(prior, function(mu, nu) c(a = mu, b = nu) + rnorm(2, sd = 0.3),
      n_simulations = n, density_estimator = "linear_gaussian")
}

test_that("a torch-free fit round-trips unchanged", {
  set.seed(1)
  fit <- fit_toy()
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  expect_identical(save_npe(fit, path), path)
  expect_true(file.exists(path))
  fit2 <- load_npe(path)

  expect_s3_class(fit2, "nsbi_npe")
  expect_equal(fit2$param_names, c("mu", "nu"))
  expect_equal(fit2$x_names, c("a", "b"))
  expect_equal(fit2$n_simulations, fit$n_simulations)
  expect_equal(fit2$std_theta, fit$std_theta)

  x_obs <- c(0.4, -0.6)
  set.seed(2); before <- sample(posterior(fit, x_obs = x_obs), 500)
  set.seed(2); after <- sample(posterior(fit2, x_obs = x_obs), 500)
  expect_equal(before, after)
})

test_that("load_npe rejects a file it did not write", {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(list(a = 1), path)
  expect_error(load_npe(path), "was not written by save_npe")
})

test_that("save_npe only accepts a fit", {
  expect_error(save_npe(list(), tempfile()), "nsbi_npe")
  expect_error(save_npe(fit_toy(100), c("a", "b")), "single file path")
})

test_that("a fit whose network died points at save_npe(), not at a torch error", {
  # stand in for what readRDS() returns: a module whose external pointer is nil
  dead <- structure(list(), class = "nsbi_dead_net")
  registerS3method("$", "nsbi_dead_net",
                   function(x, name) stop("external pointer is not valid"))
  expect_false(torch_net_alive(dead))

  fit <- fit_toy(200)
  fit$de$net <- dead
  expect_error(posterior(fit, x_obs = c(0, 0)), "save_npe")
  expect_error(save_npe(fit, tempfile()), "dangling external pointer")
  expect_output(print(fit), "network unusable")
})

test_that("de_rebuild_net refuses an estimator it cannot rebuild", {
  expect_error(de_rebuild_net(structure(list(), class = "nsbi_de_lingauss")),
               "Cannot rebuild")
})

test_that("a torch-backed fit survives save_npe() but not saveRDS()", {
  skip_if_no_torch()
  set.seed(3)
  torch::torch_manual_seed(3)
  prior <- prior_normal(mean = c(mu = 0, nu = 0), sd = 1)
  fit <- npe(prior, function(mu, nu) c(mu, nu) + rnorm(2, sd = 0.3),
             n_simulations = 400, density_estimator = "maf",
             max_epochs = 5L, hidden = c(16L, 16L), n_transforms = 2L)

  x_obs <- c(0.3, -0.4)
  theta_grid <- matrix(c(0, 0, 0.3, -0.4, 1, 1), ncol = 2, byrow = TRUE)
  before <- log_prob(posterior(fit, x_obs = x_obs), theta_grid)

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  save_npe(fit, path)
  fit2 <- load_npe(path)
  after <- log_prob(posterior(fit2, x_obs = x_obs), theta_grid)
  expect_equal(after, before, tolerance = 1e-6)

  # the failure save_npe() exists to avoid
  rds <- tempfile(fileext = ".rds")
  on.exit(unlink(rds), add = TRUE)
  saveRDS(fit, rds)
  broken <- readRDS(rds)
  expect_error(posterior(broken, x_obs = x_obs), "save_npe")
})
