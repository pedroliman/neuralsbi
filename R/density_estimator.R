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
#' @references Papamakarios, G., Pavlakou, T. and Murray, I. (2017). Masked
#'   Autoregressive Flow for Density Estimation. *NeurIPS*.
#'   \doi{10.48550/arXiv.1705.07057}
#'
#'   Durkan, C., Bekasov, A., Murray, I. and Papamakarios, G. (2019). Neural
#'   Spline Flows. *NeurIPS*. \doi{10.48550/arXiv.1906.04032}
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
#' `mu`, just on the other operand), moves both to the net's device (see
#' [net_device()]), and evaluates `log_prob_fn` under `with_no_grad()`.
#' `log_prob_fn` is the per-estimator tensor function --
#' `mdn_log_prob_tensor()`, `maf_log_prob_tensor()` or `nsf_log_prob_tensor()`.
#' The result comes back to CPU before it leaves torch, so the rest of the
#' pipeline (and everything downstream of `de_log_prob()`) stays device-agnostic
#' plain R.
#' @keywords internal
de_log_prob_torch <- function(de, theta, x, log_prob_fn) {
  theta <- as_theta_matrix(theta, de$dim_theta)
  x <- as_theta_matrix(x, de$dim_x)
  if (nrow(x) == 1L && nrow(theta) > 1L) {
    x <- matrix(x, nrow = nrow(theta), ncol = ncol(x), byrow = TRUE)
  }
  dev <- net_device(de$net)
  tt <- torch::torch_tensor(theta, dtype = torch::torch_float(), device = dev)
  xt <- torch::torch_tensor(x, dtype = torch::torch_float(), device = dev)
  torch::with_no_grad({
    lp <- log_prob_fn(de$net, tt, xt)$to(device = "cpu", dtype = torch::torch_float64())
    as.numeric(lp)
  })
}

#' Shared tensor plumbing behind `de_sample.nsbi_de_maf` and `de_sample.nsbi_de_nsf`
#'
#' Takes the single conditioning row, replicates it to `n` rows, draws a
#' standard-normal base sample, and inverts the flow -- all on the net's
#' device (see [net_device()]), then brings the draws back to CPU as a plain R
#' matrix. `inverse_fn` is the per-flow inverse -- `maf_inverse()` or
#' `nsf_inverse()`. The MDN has no inverse to share here; it samples its
#' mixture directly.
#' @keywords internal
de_sample_flow <- function(de, x, n, inverse_fn) {
  x <- as_theta_matrix(x, de$dim_x)[1, , drop = FALSE]
  xrep <- matrix(x, nrow = n, ncol = de$dim_x, byrow = TRUE)
  dev <- net_device(de$net)
  xt <- torch::torch_tensor(xrep, dtype = torch::torch_float(), device = dev)
  u <- torch::torch_randn(c(n, de$dim_theta), device = dev)
  torch::with_no_grad({
    draws <- inverse_fn(de$net, u, xt)$to(device = "cpu", dtype = torch::torch_float64())
    torch::as_array(draws)
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
#'
#' `device` is a raw keyword here (`"cpu"`, `"cuda"`, `"mps"`, `"gpu"` or
#' `"auto"`) -- resolving it to an actual, available device needs `torch`
#' loaded, so that happens inside `train_conditional_de()`, after its own
#' argument checks (`check_train_controls()`) have already run without
#' needing `torch` at all. The *resolved* string comes back on
#' `train_conditional_de()`'s return value and is stored on the returned
#' estimator (never a torch device object, which would not survive
#' `saveRDS()`) so `posterior()`/`sample()` can see what it was actually fit
#' with.
#' @inheritParams train_conditional_de
#' @keywords internal
fit_torch_de <- function(theta, x, build_net_fn, log_prob_fn, class, arch,
                         max_epochs, batch_size, lr, validation_fraction,
                         patience, n_restarts, clip_grad_norm, embedding,
                         seed, verbose, device = "cpu", min_val_rows = 1L) {
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
    seed = seed, verbose = verbose, device = device, min_val_rows = min_val_rows
  )

  structure(
    c(list(net = trained$net, dim_theta = dim_theta, dim_x = dim_x),
      arch,
      list(embedding = embedding, best_val_loss = trained$best_val_loss,
           history = trained$history, device = trained$device)),
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
  # Ridge-regularized least squares: B = (X'X + rI)^-1 X'theta. The ridge is
  # relative to each column's own scale rather than an absolute 1e-6. On
  # standardized data the two are the same thing, but this estimator also runs
  # under standardize = FALSE, and there an absolute ridge is only small if the
  # data happen to be O(1): a target column with variance 2.5e-07 had 1e-06
  # added to it and came back five times too wide.
  XtX <- crossprod(X)
  diag(XtX) <- diag(XtX) + ridge * ridge_scale(diag(XtX))
  B <- solve(XtX, crossprod(X, theta))   # (q+1) x p
  mu <- X %*% B
  resid <- theta - mu
  Sigma <- crossprod(resid) / max(n - ncol(X), 1)
  # The target's own variance, not the residual's, is what sets the scale here.
  # A noiseless model leaves no residual at all, and that is exactly when the
  # ridge has to supply a positive diagonal for chol() to work with.
  diag(Sigma) <- diag(Sigma) +
    ridge * ridge_scale(apply(theta, 2, stats::var))
  verbose_cat(verbose, sprintf(
    "[linear_gaussian] fitted on %d sims, %d params, %d data dims\n",
    n, p, ncol(x)))
  structure(
    list(B = B, Sigma = Sigma, chol = chol(Sigma), dim_theta = p,
         dim_x = ncol(x)),
    class = c("nsbi_de_lingauss", "nsbi_de")
  )
}

#' Per-column scale a relative ridge is measured against
#'
#' One scale per column, so a matrix whose columns span several orders of
#' magnitude is regularized by the same *relative* amount everywhere. A single
#' shared scale -- the mean or the maximum of the diagonal -- is no better than
#' an absolute ridge here: it is set by the largest column and swamps the
#' smallest.
#'
#' A column with no scale of its own borrows the largest one that has any, and
#' if nothing does the scale is 1 and the ridge acts absolutely. That is the
#' degenerate case the ridge exists for, and it is the only one where the
#' answer is arbitrary.
#'
#' @param scale Per-column scales, in the units of the matrix diagonal.
#' @keywords internal
ridge_scale <- function(scale) {
  s <- as.numeric(scale)
  bad <- !is.finite(s) | s <= 0
  if (all(bad)) return(rep(1, length(s)))
  s[bad] <- max(s[!bad])
  s
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
