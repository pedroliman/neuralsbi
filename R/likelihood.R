#' Evaluate a surrogate likelihood
#'
#' `log_lik()` evaluates the likelihood learned by [nle()]:
#' \eqn{\log q_\phi(x \mid \theta)}, in the original units of both.
#'
#' Rows of `x` are treated as **independent observations from the same
#' parameter**, so by default the result sums over them,
#' \eqn{\sum_i \log q_\phi(x_i \mid \theta)}. That sum is the whole point of
#' NLE: the estimator is trained on one observation at a time, and the number
#' of observations you condition on afterwards is free.
#'
#' @param fit An `nsbi_nle` fit from [nle()].
#' @param theta Parameter values: a numeric vector (one parameter set) or an
#'   `n_theta x dim_theta` matrix.
#' @param x Observed data: a numeric vector (one observation) or an
#'   `n_obs x dim_x` matrix whose rows are independent observations.
#' @param sum_iid Sum the log-density over the rows of `x` (the default). Set
#'   `FALSE` to get the per-observation values instead.
#' @param max_batch Largest number of `(theta, x)` pairs evaluated in one call
#'   to the estimator. Only affects memory and speed.
#' @param ... Unused, for S3 consistency.
#'
#' @return With `sum_iid = TRUE`, a numeric vector with one entry per row of
#'   `theta`. With `sum_iid = FALSE`, an `n_theta x n_obs` matrix.
#'
#' @seealso [likelihood_fn()] for a closure over a fixed observation,
#'   [posterior()] to turn the likelihood into posterior draws.
#'
#' @examples
#' prior <- prior_uniform(c(mu = -3), c(mu = 3))
#' fit <- nle(prior, function(mu) c(y = rnorm(1, mu, 0.5)),
#'            n_simulations = 1000, density_estimator = "linear_gaussian")
#'
#' x_obs <- matrix(rnorm(20, mean = 1, sd = 0.5), ncol = 1)
#' grid <- matrix(seq(-2, 2, length.out = 5), ncol = 1)
#' log_lik(fit, grid, x_obs)
#' @export
log_lik <- function(fit, theta, x, ...) UseMethod("log_lik")

#' @rdname log_lik
#' @export
log_lik.nsbi_nle <- function(fit, theta, x, sum_iid = TRUE,
                             max_batch = 1e5, ...) {
  check_fit_alive(fit)
  # Both arguments are matrices with a required width, so a bare "expected 2
  # columns" would leave the user guessing which one is wrong.
  theta <- as_lik_matrix(theta, fit$dim_theta, "theta", "one parameter per column")
  x <- as_lik_matrix(x, fit$dim_x, "x",
                     "one row per independent observation")

  theta_z <- apply_standardizer(fit$std_theta, theta)
  x_z <- apply_standardizer(fit$std_x, x)
  n_theta <- nrow(theta_z)
  n_obs <- nrow(x_z)

  out <- de_log_lik_iid(fit$de, x_z, theta_z, max_batch = max_batch)

  # de_log_prob works in standardized x space; this puts it back in the units
  # the simulator returned.
  out <- out + standardizer_log_jac(fit$std_x)

  if (!isTRUE(sum_iid)) {
    return(out)
  }
  stats::setNames(rowSums(out), rownames(theta))
}

#' Log-density of many observations under many parameter values
#'
#' The cross product `n_theta x n_obs`, in standardized space. This is the hot
#' loop behind [log_lik()] and behind every MCMC step, so it gets its own
#' generic rather than going through [de_log_prob()] for every pair.
#'
#' The default expands the cross product and makes one batched call, which is
#' what a flow needs: its transforms depend on the observation as well as the
#' parameter, so there is nothing to reuse between observations. Estimators
#' whose conditional distribution depends on the parameter *alone* -- the MDN
#' and the linear-Gaussian baseline -- override this and compute that
#' distribution once per parameter, which turns the i.i.d. sum from
#' `n_theta * n_obs` network passes into `n_theta` of them. With a few thousand
#' observations that is the difference between usable and not.
#'
#' @param de A fitted density estimator.
#' @param x Standardized observations, `n_obs x dim_x`.
#' @param theta Standardized parameters, `n_theta x dim_theta`.
#' @param max_batch Largest number of pairs evaluated at once.
#' @return An `n_theta x n_obs` matrix of log-densities.
#' @keywords internal
de_log_lik_iid <- function(de, x, theta, max_batch = 1e5) {
  UseMethod("de_log_lik_iid")
}

#' @export
de_log_lik_iid.default <- function(de, x, theta, max_batch = 1e5) {
  n_theta <- nrow(theta)
  n_obs <- nrow(x)
  per_block <- max(1L, floor(max_batch / max(n_obs, 1L)))
  out <- matrix(0, nrow = n_theta, ncol = n_obs)
  for (s in seq.int(1L, n_theta, by = per_block)) {
    idx <- seq.int(s, min(s + per_block - 1L, n_theta))
    cond <- theta[rep(idx, each = n_obs), , drop = FALSE]
    target <- x[rep(seq_len(n_obs), times = length(idx)), , drop = FALSE]
    out[idx, ] <- matrix(de_log_prob(de, target, cond),
                         nrow = length(idx), ncol = n_obs, byrow = TRUE)
  }
  out
}

#' @export
de_log_lik_iid.nsbi_de_lingauss <- function(de, x, theta, max_batch = 1e5) {
  # One conditional mean per parameter, then every observation scored against
  # it. No loop over observations at all.
  mu <- lingauss_mean(de, theta)
  lp <- vapply(seq_len(nrow(mu)),
               function(i) dmvnorm_chol(x, mu[i, ], de$chol, log = TRUE),
               numeric(nrow(x)))
  # vapply returns a vector, not a one-row matrix, when there is a single
  # observation, so the shape is set explicitly rather than by transposing.
  matrix(lp, nrow = nrow(mu), ncol = nrow(x), byrow = TRUE)
}

#' @export
de_log_lik_iid.nsbi_de_mdn <- function(de, x, theta, max_batch = 1e5) {
  # The MLP maps theta to the mixture parameters and never sees x, so the
  # forward pass runs once per parameter and every observation is scored
  # against the resulting mixture.
  n_theta <- nrow(theta)
  n_obs <- nrow(x)
  K <- de$n_components
  p <- de$dim_theta                    # the estimator's target dimension

  tt <- torch::torch_tensor(theta, dtype = torch::torch_float())
  xt <- torch::torch_tensor(x, dtype = torch::torch_float())
  torch::with_no_grad({
    params <- de$net(tt)
    L <- mdn_build_tril(de$net, params$tril_flat)                  # (T,K,p,p)
    log_w <- torch::nnf_log_softmax(params$logits, dim = 2)        # (T,K)
    diag_L <- torch::torch_diagonal(L, dim1 = 3, dim2 = 4)
    logdet <- 2 * torch::torch_log(diag_L)$sum(dim = 3)            # (T,K)

    # Chunk the observations so the (T, K, p, n) intermediate stays bounded.
    per_chunk <- max(1L, floor(max_batch / max(n_theta * K * p, 1L)))
    pieces <- list()
    for (s in seq.int(1L, n_obs, by = per_chunk)) {
      idx <- seq.int(s, min(s + per_chunk - 1L, n_obs))
      xs <- xt[idx, , drop = FALSE]                                # (n,p)
      # (T,K,n,p): every observation minus every component mean
      diff <- xs$unsqueeze(1)$unsqueeze(1) - params$means$unsqueeze(3)
      z <- torch::linalg_solve_triangular(L, diff$transpose(3, 4), upper = FALSE)
      quad <- z$pow(2)$sum(dim = 3)                                # (T,K,n)
      comp <- -0.5 * (p * log(2 * pi) + logdet$unsqueeze(3) + quad)
      lp <- torch::torch_logsumexp(log_w$unsqueeze(3) + comp, dim = 2)
      pieces[[length(pieces) + 1L]] <-
        torch::as_array(lp$to(dtype = torch::torch_float64()))
    }
    matrix(unlist(pieces), nrow = n_theta)
  })
}

#' Coerce one of `log_lik()`'s two matrix arguments, naming it if it is wrong
#' @keywords internal
as_lik_matrix <- function(value, d, arg, what) {
  tryCatch(
    as_theta_matrix(value, d),
    error = function(e) {
      stop(sprintf("`%s` must have %d columns (%s).", arg, d, what),
           call. = FALSE)
    }
  )
}

#' A surrogate likelihood as a plain R function
#'
#' Fixes an observation and returns `function(theta)`, so the learned
#' likelihood can be handed to anything in R that wants a log-density:
#' [stats::optim()], an MCMC package, an importance sampler, a profile
#' likelihood. The returned function is vectorized -- pass a matrix of
#' parameter values and get one log-likelihood per row.
#'
#' @param fit An `nsbi_nle` fit from [nle()].
#' @param x_obs The observation to condition on. Rows are independent
#'   observations; the log-likelihood sums over them.
#' @param ... Passed to [log_lik()].
#'
#' @return A function of `theta` returning a numeric vector of
#'   log-likelihoods, with the observation attached as attribute `x_obs`.
#'
#' @examples
#' prior <- prior_uniform(c(mu = -3), c(mu = 3))
#' fit <- nle(prior, function(mu) c(y = rnorm(1, mu, 0.5)),
#'            n_simulations = 1000, density_estimator = "linear_gaussian")
#'
#' loglik <- likelihood_fn(fit, matrix(rnorm(20, 1, 0.5), ncol = 1))
#' optimize(function(m) loglik(m), c(-3, 3), maximum = TRUE)$maximum
#' @export
likelihood_fn <- function(fit, x_obs, ...) {
  stopifnot(inherits(fit, "nsbi_nle"))
  x_obs <- as_theta_matrix(x_obs, fit$dim_x)
  force(x_obs)
  dots <- list(...)
  f <- function(theta) {
    do.call(log_lik, c(list(fit = fit, theta = theta, x = x_obs), dots))
  }
  attr(f, "x_obs") <- x_obs
  f
}

#' Unnormalized log posterior of an NLE fit
#'
#' \eqn{\log q_\phi(x \mid \theta) + \log p(\theta)}, returning `-Inf` outside
#' the prior support. This is the potential the MCMC samplers target.
#' @keywords internal
nle_potential <- function(fit, x_obs, max_batch = 1e5) {
  prior <- fit$prior
  # prior_custom() without a log_prob_fn returns NA rather than nothing, so the
  # only way to find out is to ask it. Better here than as a puzzling
  # initialization failure a few hundred lines later.
  probe <- tryCatch(prior$log_prob(sample_prior(prior, 2L)),
                    error = function(e) NA_real_)
  if (is.null(prior$log_prob) || all(is.na(probe))) {
    stop("Sampling an NLE posterior needs a prior log-density, and this prior ",
         "does not have one.\nRebuild it with prior_custom(..., log_prob_fn = ).",
         call. = FALSE)
  }
  x_obs <- as_theta_matrix(x_obs, fit$dim_x)
  bounded <- !is.null(prior$lower) || !is.null(prior$upper)
  function(theta) {
    theta <- as_theta_matrix(theta, fit$dim_theta)
    lp <- as.numeric(prior$log_prob(theta))
    if (bounded) lp[!within_support(prior, theta)] <- -Inf
    ok <- is.finite(lp)
    if (!any(ok)) return(lp)
    # Only evaluate the network where the prior gives the point any mass; that
    # is a large saving once a chain starts probing a bounded edge.
    lp[ok] <- lp[ok] + log_lik(fit, theta[ok, , drop = FALSE], x_obs,
                               max_batch = max_batch)
    lp
  }
}
