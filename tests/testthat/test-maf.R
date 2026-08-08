# Masked Autoregressive Flow: masks, invertibility, and posterior accuracy.

test_that("MADE masks enforce the autoregressive property", {
  m <- made_masks(dim_theta = 3, dim_x = 2, hidden = c(8, 8))
  # x inputs (degree 0) reach every hidden unit
  expect_true(all(m$hidden[[1]][, 4:5] == 1))
  # output d must not see hidden units of degree >= d: output 1 sees nothing
  # from theta (all hidden degrees >= 1)
  expect_true(all(m$out[1, ] == 0))
  # for p = 1 the (only) output depends on x alone
  m1 <- made_masks(dim_theta = 1, dim_x = 2, hidden = c(4))
  expect_true(all(m1$hidden[[1]][, 1] == 0))  # theta input cut off
  expect_true(all(m1$out == 1))               # x-driven hidden units visible
})

test_that("MAF forward and inverse are consistent (round trip)", {
  skip_if_no_torch()
  torch::torch_manual_seed(5)
  net <- maf_module(dim_x = 2, dim_theta = 3, n_transforms = 3,
                    hidden = c(16, 16))()
  # give the (identity-initialized) output heads nonzero weights
  for (p in net$parameters) torch::with_no_grad(p$add_(torch::torch_randn_like(p) * 0.1))
  theta <- torch::torch_randn(c(20, 3))
  x <- torch::torch_randn(c(20, 2))
  torch::with_no_grad({
    fw <- maf_forward(net, theta, x)
    back <- maf_inverse(net, fw$u, x)
  })
  expect_lt(max(abs(torch::as_array(back - theta))), 1e-4)
})

test_that("MAF log_prob integrates the change of variables correctly", {
  skip_if_no_torch()
  # with identity-initialized transforms, log q(theta|x) is the standard normal
  net <- maf_module(dim_x = 1, dim_theta = 2, n_transforms = 2,
                    hidden = c(8))()
  theta <- matrix(c(0, 0, 1, -1), nrow = 2, byrow = TRUE)
  x <- matrix(0, nrow = 2, ncol = 1)
  lp <- de_log_prob(structure(list(net = net, dim_theta = 2, dim_x = 1),
                              class = c("nsbi_de_maf", "nsbi_de")),
                    theta, x)
  expect_equal(lp, rowSums(dnorm(theta, log = TRUE)), tolerance = 1e-5)
})

test_that("MAF posterior is close to the analytic linear-Gaussian posterior", {
  skip_if_no_torch()
  set.seed(6)
  d <- 2; sigma <- 0.5
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + rnorm(length(theta), sd = sigma)
  fit <- npe(prior, simulator, n_simulations = 8000, density_estimator = "maf",
             n_transforms = 3L, hidden = c(50L, 50L), max_epochs = 300L,
             seed = 6)
  x_obs <- c(1.0, -0.5)
  post <- posterior(fit, x_obs = x_obs)
  draws <- sample(post, 10000)

  truth <- analytic_gauss_posterior(x_obs, sigma, d)
  mu <- truth$mu; Sigma <- truth$Sigma
  expect_equal(colMeans(draws), mu, tolerance = 0.1)
  expect_equal(apply(draws, 2, sd), sqrt(diag(Sigma)), tolerance = 0.1)
})
