#' Conditional density estimators
#'
#' A conditional density estimator learns \eqn{q_\phi(\theta \mid x)}. In
#' `neuralsbi` every estimator is trained in *standardized* space and exposes two
#' generics:
#'
#' * `de_log_prob(de, theta, x)` -- log density of `theta` given `x`
#' * `de_sample(de, x, n)` -- draw `n` parameter vectors given a single `x`
#'
#' The contract is really `q(target | condition)`: it makes no assumption
#' about which of the two arguments is the parameter. [npe()] calls it with
#' `theta` as the target and `x` as the condition, learning the posterior.
#' [nle()] swaps the two, learning the likelihood instead, with the same
#' estimators and the same two generics.
#'
#' Four estimators ship today:
#'
#' * `"maf"` -- a Masked Autoregressive Flow (Papamakarios et al., 2017), a
#'   stack of invertible autoregressive transforms with an exact
#'   change-of-variables density. This is the default, matching Python
#'   `sbi`, and requires the `torch` back end.
#' * `"nsf"` -- a Neural Spline Flow (Durkan et al., 2019): the same
#'   autoregressive structure as the MAF, but with a monotonic
#'   rational-quadratic spline transform in place of MAF's affine one, which
#'   handles sharply non-Gaussian posteriors better. Requires `torch`.
#' * `"mdn"` -- a Mixture Density Network (neural network -> Gaussian
#'   mixture). Requires `torch`.
#' * `"linear_gaussian"` -- a closed-form conditional Gaussian baseline
#'   (least-squares mean, residual covariance). No neural network, no `torch`.
#'   It is exact for linear-Gaussian simulators and doubles as a fast baseline
#'   and a regression-test oracle.
#'
#' @name density_estimator
NULL

#' @keywords internal
de_log_prob <- function(de, theta, x) UseMethod("de_log_prob")

#' @keywords internal
de_sample <- function(de, x, n) UseMethod("de_sample")

# ---- shared tensor plumbing for the torch estimators (MDN, MAF, NSF) ------

#' Shared tensor plumbing behind every neural `de_log_prob.*` method
#'
#' Coerces `theta` and `x` to matrices, broadcasts a single-row `x` up to
#' `theta`'s row count (the same broadcast [lingauss_mean()]'s caller does on
#' `mu`, just on the other operand), moves both to torch, and evaluates
#' `log_prob_fn` under `with_no_grad()`. `log_prob_fn` is the per-estimator
#' tensor function -- `mdn_log_prob_tensor()`, `maf_log_prob_tensor()` or
#' `nsf_log_prob_tensor()`.
#' @keywords internal
de_log_prob_torch <- function(de, theta, x, log_prob_fn) {
  theta <- as_theta_matrix(theta, de$dim_theta)
  x <- as_theta_matrix(x, de$dim_x)
  if (nrow(x) == 1L && nrow(theta) > 1L) {
    x <- matrix(x, nrow = nrow(theta), ncol = ncol(x), byrow = TRUE)
  }
  tt <- torch::torch_tensor(theta, dtype = torch::torch_float())
  xt <- torch::torch_tensor(x, dtype = torch::torch_float())
  torch::with_no_grad({
    as.numeric(log_prob_fn(de$net, tt, xt)$to(dtype = torch::torch_float64()))
  })
}

#' Shared tensor plumbing behind `de_sample.nsbi_de_maf` and `de_sample.nsbi_de_nsf`
#'
#' Takes the single conditioning row, replicates it to `n` rows, draws a
#' standard-normal base sample, and inverts the flow. `inverse_fn` is the
#' per-flow inverse -- `maf_inverse()` or `nsf_inverse()`. The MDN has no
#' inverse to share here; it samples its mixture directly.
#' @keywords internal
de_sample_flow <- function(de, x, n, inverse_fn) {
  x <- as_theta_matrix(x, de$dim_x)[1, , drop = FALSE]
  xrep <- matrix(x, nrow = n, ncol = de$dim_x, byrow = TRUE)
  xt <- torch::torch_tensor(xrep, dtype = torch::torch_float())
  u <- torch::torch_randn(c(n, de$dim_theta))
  torch::with_no_grad({
    torch::as_array(inverse_fn(de$net, u, xt)$to(dtype = torch::torch_float64()))
  })
}

#' Shared body behind `fit_mdn()`, `fit_maf()` and `fit_nsf()`
#'
#' Coerces `theta` and `x`, builds the net from `build_net_fn(dim_x,
#' dim_theta)` now that both are known, trains it with
#' [train_conditional_de()], and packages the result into a fitted `nsbi_de`
#' object.
#'
#' `arch` carries the architecture fields specific to the caller --
#' `n_components`/`hidden` for the MDN, `n_transforms`/`hidden` for the MAF,
#' `n_transforms`/`hidden`/`n_bins`/`tail_bound` for the NSF -- and is spliced
#' into the returned list ahead of `embedding`, matching the field order each
#' estimator returned before this helper existed. This helper never needs to
#' know what `arch`'s fields are.
#' @keywords internal
fit_torch_de <- function(theta, x, build_net_fn, log_prob_fn, class, arch,
                         max_epochs, batch_size, lr, validation_fraction,
                         patience, n_restarts, clip_grad_norm, embedding,
                         seed, verbose) {
  theta <- as_theta_matrix(theta)
  x <- as_theta_matrix(x)
  dim_theta <- ncol(theta)
  dim_x <- ncol(x)

  trained <- train_conditional_de(
    build_net = function() build_net_fn(dim_x, dim_theta),
    log_prob_fn = log_prob_fn,
    theta = theta, x = x,
    max_epochs = max_epochs, batch_size = batch_size, lr = lr,
    validation_fraction = validation_fraction, patience = patience,
    n_restarts = n_restarts, clip_grad_norm = clip_grad_norm,
    seed = seed, verbose = verbose
  )

  structure(
    c(list(net = trained$net, dim_theta = dim_theta, dim_x = dim_x),
      arch,
      list(embedding = embedding, best_val_loss = trained$best_val_loss,
           history = trained$history)),
    class = c(class, "nsbi_de")
  )
}

# ---- linear-Gaussian conditional estimator (pure R) -----------------------

#' @keywords internal
fit_linear_gaussian <- function(theta, x, ridge = 1e-6, verbose = FALSE) {
  theta <- as_theta_matrix(theta)
  x <- as_theta_matrix(x)
  n <- nrow(theta)
  p <- ncol(theta)
  X <- cbind(1, x)                       # design matrix with intercept
  # Ridge-regularized least squares: B = (X'X + rI)^-1 X'theta
  XtX <- crossprod(X)
  diag(XtX) <- diag(XtX) + ridge
  B <- solve(XtX, crossprod(X, theta))   # (q+1) x p
  mu <- X %*% B
  resid <- theta - mu
  Sigma <- crossprod(resid) / max(n - ncol(X), 1)
  diag(Sigma) <- diag(Sigma) + ridge
  verbose_cat(verbose, sprintf(
    "[linear_gaussian] fitted on %d sims, %d params, %d data dims\n",
    n, p, ncol(x)))
  structure(
    list(B = B, Sigma = Sigma, chol = chol(Sigma), dim_theta = p,
         dim_x = ncol(x)),
    class = c("nsbi_de_lingauss", "nsbi_de")
  )
}

#' Conditional mean of the linear-Gaussian estimator
#'
#' `lingauss_mean()` is the only place `x` meets `de$B`, so it is also the
#' place the width of `x` is checked. It passes `de$dim_x` to
#' [as_theta_matrix()] for the reason every neural estimator passes its own: a
#' wrong-width `x` otherwise gets as far as the matrix product and is reported
#' as "non-conformable arguments", which names neither the argument nor the
#' width expected of it. An estimator fitted before `dim_x` was recorded has
#' `NULL` here and keeps the old unchecked behaviour.
#'
#' @param de A fitted `nsbi_de_lingauss` object.
#' @param x The conditioning variable: a numeric vector, matrix or data frame
#'   with `de$dim_x` columns.
#' @return An `nrow(x) x de$dim_theta` matrix of conditional means.
#' @keywords internal
lingauss_mean <- function(de, x) {
  x <- as_theta_matrix(x, de$dim_x)
  cbind(1, x) %*% de$B
}

#' @export
de_log_prob.nsbi_de_lingauss <- function(de, theta, x) {
  theta <- as_theta_matrix(theta, de$dim_theta)
  mu <- lingauss_mean(de, x)
  if (nrow(mu) == 1L && nrow(theta) > 1L) {
    mu <- matrix(mu, nrow = nrow(theta), ncol = ncol(mu), byrow = TRUE)
  }
  dmvnorm_chol(theta, mu, de$chol)
}

#' @export
de_sample.nsbi_de_lingauss <- function(de, x, n) {
  mu <- lingauss_mean(de, x)[1, ]
  z <- matrix(stats::rnorm(n * de$dim_theta), nrow = n)
  sweep(z %*% de$chol, 2, mu, `+`)
}

# ---- small multivariate-normal helpers ------------------------------------

#' Multivariate normal log density using a precomputed upper-Cholesky factor
#' (`R` such that `Sigma = t(R) %*% R`, i.e. `chol(Sigma)`).
#'
#' Always returns the log density: every call site wants `log_prob()`'s
#' contract, none of `dnorm()`'s `log = FALSE`, so there is no `log` argument
#' to forget to set.
#' @keywords internal
dmvnorm_chol <- function(x, mean, R) {
  x <- as_theta_matrix(x)
  if (is.null(dim(mean))) mean <- matrix(mean, nrow = nrow(x), ncol = ncol(x),
                                         byrow = TRUE)
  d <- ncol(x)
  dev <- x - mean
  # Solve R' z = dev'  =>  quadratic form = colSums(z^2)
  z <- backsolve(R, t(dev), transpose = TRUE)
  quad <- colSums(z^2)
  logdet <- 2 * sum(log(diag(R)))
  -0.5 * (d * log(2 * pi) + logdet + quad)
}
