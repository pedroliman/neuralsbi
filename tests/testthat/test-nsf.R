# Neural Spline Flow: spline math, invertibility, and posterior accuracy.

test_that("rational-quadratic spline round-trips and is identity outside", {
  skip_if_no_torch()
  torch::torch_manual_seed(12)
  N <- 200L; K <- 8L
  xin <- torch::torch_cat(list(torch::torch_rand(N) * 5 - 2.5,  # inside
                               torch::torch_tensor(c(-5, 5, 3.5))))
  w <- torch::torch_randn(c(xin$shape[1], K))
  h <- torch::torch_randn(c(xin$shape[1], K))
  d <- torch::torch_randn(c(xin$shape[1], K - 1L))
  fw <- rq_spline(xin, w, h, d, inverse = FALSE, tail_bound = 3)
  bk <- rq_spline(fw$outputs, w, h, d, inverse = TRUE, tail_bound = 3)
  expect_lt(max(abs(torch::as_array(bk$outputs - xin))), 1e-4)
  # forward and inverse log-dets cancel
  expect_lt(max(abs(torch::as_array(fw$logdet + bk$logdet))), 1e-4)
  # identity with zero log-det outside the tail bound
  outside <- torch::as_array(xin$abs() > 3)
  expect_equal(torch::as_array(fw$outputs)[outside],
               torch::as_array(xin)[outside], tolerance = 1e-6)
  expect_equal(torch::as_array(fw$logdet)[outside],
               rep(0, sum(outside)), tolerance = 1e-6)
})

test_that("NSF flow forward/inverse round trip through the full stack", {
  skip_if_no_torch()
  torch::torch_manual_seed(13)
  p_dim <- 3L
  net <- nsf_module(dim_x = 2, dim_theta = p_dim, n_transforms = 3,
                    hidden = c(16), n_bins = 5, tail_bound = 3)()
  for (p in net$parameters) torch::with_no_grad(p$add_(torch::torch_randn_like(p) * 0.05))
  theta <- torch::torch_randn(c(15, p_dim))
  x <- torch::torch_randn(c(15, 2))
  torch::with_no_grad({
    lp <- nsf_log_prob_tensor(net, theta, x)
    # forward through the stack (mirrors nsf_log_prob_tensor), then invert
    z <- theta
    rev_idx <- torch::torch_tensor(rev(seq_len(p_dim)), dtype = torch::torch_long())
    for (k in seq_len(net$n_transforms)) {
      if (k > 1L) z <- z[, rev_idx, drop = FALSE]
      z <- nsf_apply(net, net$mades[[k]], z, x, inverse = FALSE)$z
    }
    back <- nsf_inverse(net, z, x)
  })
  expect_true(all(is.finite(torch::as_array(lp))))
  expect_lt(max(abs(torch::as_array(back - theta))), 1e-4)
})

test_that("NSF trains end to end and yields a sane linear-Gaussian posterior", {
  # A lightweight smoke test: the flow's correctness is covered rigorously by
  # the round-trip / tail-identity / log-det tests above; here we only confirm
  # that a trained NSF produces a posterior in the right place with finite
  # densities. Tight analytic parity is checked (non-CI) via inst/benchmarks/.
  skip_if_no_torch()
  set.seed(14)
  d <- 2; sigma <- 0.5
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  simulator <- function(theta) theta + matrix(rnorm(length(theta), sd = sigma),
                                              nrow = nrow(theta))
  fit <- npe(prior, simulator, n_simulations = 1200, density_estimator = "nsf",
             n_transforms = 2L, hidden = c(24L, 24L), max_epochs = 50L,
             seed = 14)
  x_obs <- c(1.0, -0.5)
  post <- posterior(fit, x_obs = x_obs)
  draws <- sample(post, 5000)

  prec <- diag(d) + diag(d) / sigma^2
  Sigma <- solve(prec)
  mu <- as.numeric(Sigma %*% (x_obs / sigma^2))
  # posterior mean in the right neighborhood; sd in a plausible band
  expect_lt(max(abs(colMeans(draws) - mu)), 0.2)
  expect_true(all(apply(draws, 2, sd) > 0.3 & apply(draws, 2, sd) < 0.7))
  # densities are finite
  lp <- log_prob(post, draws[1:20, ], normalize = FALSE)
  expect_true(all(is.finite(lp)))
})
