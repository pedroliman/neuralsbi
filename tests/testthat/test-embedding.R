# Embedding networks: spec construction, wiring into the estimators, and that a
# jointly trained summary net still recovers the analytic posterior.

test_that("embedding_mlp builds a torch-free spec", {
  emb <- embedding_mlp(output_dim = 8, hidden = c(32, 16))
  expect_s3_class(emb, "nsbi_embedding")
  expect_identical(emb$output_dim, 8L)
  expect_identical(emb$hidden, c(32L, 16L))
  # effective conditioning dim is output_dim with an embedding, dim_x without
  expect_identical(embedding_output_dim(emb, dim_x = 20), 8L)
  expect_identical(embedding_output_dim(NULL, dim_x = 20), 20)
})

test_that("embedding_mlp validates its arguments", {
  expect_error(embedding_mlp(output_dim = c(4, 5)), "output_dim")
  expect_error(embedding_mlp(output_dim = 0), "output_dim")
  expect_error(embedding_mlp(hidden = c(8, -1)), "hidden")
})

test_that("npe warns and ignores an embedding for linear_gaussian", {
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + matrix(rnorm(length(theta), sd = 0.3),
                                              nrow = nrow(theta))
  expect_warning(
    fit <- npe(prior, simulator, n_simulations = 200,
               density_estimator = "linear_gaussian",
               embedding_net = embedding_mlp(output_dim = 2)),
    "linear_gaussian"
  )
  expect_null(fit$de$embedding)
})

test_that("npe rejects a non-spec embedding_net", {
  prior <- prior_normal(mean = 0, sd = 1)
  simulator <- function(theta) theta
  expect_error(
    npe(prior, simulator, n_simulations = 10, density_estimator = "mdn",
        embedding_net = list(output_dim = 3)),
    "embedding_mlp"
  )
})

test_that("the embedding module maps raw x to output_dim features", {
  skip_if_no_torch()
  torch::torch_manual_seed(1)
  mod <- build_embedding_module(embedding_mlp(output_dim = 4, hidden = c(8)),
                                dim_x = 10)
  h <- mod(torch::torch_randn(c(7, 10)))
  expect_identical(dim(torch::as_array(h)), c(7L, 4L))
  # NULL spec is the identity embedding (no module)
  expect_null(build_embedding_module(NULL, dim_x = 10))
})

test_that("embedded estimators expose the summary net and keep raw dim_x", {
  skip_if_no_torch()
  set.seed(2)
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + matrix(rnorm(length(theta), sd = 0.4),
                                              nrow = nrow(theta))
  emb <- embedding_mlp(output_dim = 3, hidden = c(16))
  for (de in c("mdn", "maf", "nsf")) {
    fit <- npe(prior, simulator, n_simulations = 300, density_estimator = de,
               embedding_net = emb, max_epochs = 5L, seed = 2)
    expect_equal(fit$dim_x, 2)          # de_* still take the raw data
    expect_true(isTRUE(fit$de$net$has_embedding))
    # the embedding's parameters are part of the trained network
    pnames <- names(fit$de$net$state_dict())
    expect_true(any(grepl("embedding", pnames)))
    post <- posterior(fit, x_obs = c(0.5, -0.5))
    draws <- sample(post, 50)
    expect_equal(dim(draws), c(50L, 2L))
    expect_true(all(is.finite(draws)))
  }
})

test_that("a jointly trained embedding still recovers the linear-Gaussian posterior", {
  skip_if_no_torch()
  set.seed(7)
  d <- 2; sigma <- 0.5
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + matrix(rnorm(length(theta), sd = sigma),
                                              nrow = nrow(theta))
  fit <- npe(prior, simulator, n_simulations = 8000, density_estimator = "maf",
             n_transforms = 3L, hidden = c(50L, 50L),
             embedding_net = embedding_mlp(output_dim = 4, hidden = c(32)),
             max_epochs = 300L, seed = 7)
  x_obs <- c(1.0, -0.5)
  post <- posterior(fit, x_obs = x_obs)
  draws <- sample(post, 10000)

  prec <- diag(d) + diag(d) / sigma^2
  Sigma <- solve(prec)
  mu <- as.numeric(Sigma %*% (x_obs / sigma^2))
  expect_equal(colMeans(draws), mu, tolerance = 0.12)
})
