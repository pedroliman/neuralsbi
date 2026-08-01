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

  # de_log_prob works in standardized x space; the Jacobian puts it back in the
  # units the simulator returned. It is constant, so the sum over observations
  # picks it up once per observation.
  jac <- standardizer_log_jac(fit$std_x)

  if (!isTRUE(sum_iid)) {
    return(de_log_lik_iid(fit$de, x_z, theta_z, max_batch = max_batch) + jac)
  }
  ev <- de_iid_evaluator(fit$de, x_z, max_batch = max_batch)
  stats::setNames(ev(theta_z) + nrow(x_z) * jac, rownames(theta))
}

#' Log-density of many observations under many parameter values
#'
#' The cross product `n_theta x n_obs`, in standardized space. Only
#' `sum_iid = FALSE` wants that matrix; everything else -- [log_lik()]'s default
#' and every MCMC step -- wants its row sums, and gets them from
#' [de_iid_evaluator()] without ever building it.
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
  n_obs <- nrow(x)
  out <- matrix(0, nrow = nrow(theta), ncol = n_obs)
  cross_iid(de, x, theta, max_batch, function(idx, lp) {
    out[idx, ] <<- matrix(lp, nrow = length(idx), ncol = n_obs, byrow = TRUE)
  })
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
  xt <- torch::torch_tensor(x, dtype = torch::torch_float())
  pieces <- list()
  mdn_iid_blocks(de, xt, theta, max_batch, function(idx, lp) {
    pieces[[length(pieces) + 1L]] <<-
      torch::as_array(lp$to(dtype = torch::torch_float64()))
  })
  matrix(unlist(pieces), nrow = nrow(theta))
}

#' A summed i.i.d. log-likelihood with the observation held fixed
#'
#' [log_lik()] and every MCMC step ask the same question over and over: the
#' summed log-density of one fixed set of observations under a `theta` that
#' changes. `de_iid_evaluator()` returns a closure over the observations, so
#' whatever an estimator can settle once settles when the closure is built
#' rather than on every call. For the MDN that is coercing the observations to
#' a tensor, which at a few thousand rows is not a rounding error next to the
#' forward pass.
#'
#' Reducing inside the closure matters as much as the hoisting. The
#' `n_theta x n_obs` matrix is the largest object in the loop and none of it is
#' wanted, so the sum happens where the log-densities are produced and only
#' `n_theta` numbers ever cross back into R.
#'
#' Estimators need a method here only if they can beat the default, which is
#' [de_log_lik_iid()] with its row sums taken block by block.
#'
#' @inheritParams de_log_lik_iid
#' @return `function(theta)` giving one summed log-density per row of `theta`,
#'   in standardized space.
#' @keywords internal
de_iid_evaluator <- function(de, x, max_batch = 1e5) UseMethod("de_iid_evaluator")

#' @export
de_iid_evaluator.default <- function(de, x, max_batch = 1e5) {
  n_obs <- nrow(x)
  force(max_batch)
  function(theta) {
    out <- numeric(nrow(theta))
    cross_iid(de, x, theta, max_batch, function(idx, lp) {
      # `lp` runs theta-major, so a matrix of it has one column per theta.
      out[idx] <<- colSums(matrix(lp, nrow = n_obs))
    })
    out
  }
}

#' @export
de_iid_evaluator.nsbi_de_lingauss <- function(de, x, max_batch = 1e5) {
  function(theta) {
    mu <- lingauss_mean(de, theta)
    vapply(seq_len(nrow(mu)),
           function(i) sum(dmvnorm_chol(x, mu[i, ], de$chol, log = TRUE)),
           numeric(1))
  }
}

#' @export
de_iid_evaluator.nsbi_de_mdn <- function(de, x, max_batch = 1e5) {
  # The observations become a tensor once, not once per MCMC step.
  xt <- torch::torch_tensor(x, dtype = torch::torch_float())
  force(max_batch)
  eager <- function(theta) {
    total <- numeric(nrow(theta))
    mdn_iid_blocks(de, xt, theta, max_batch, function(idx, lp) {
      total <<- total +
        as.numeric(lp$to(dtype = torch::torch_float64())$sum(dim = 2))
    })
    total
  }
  traced <- mdn_trace_cache(de, xt, max_batch, eager)
  if (is.null(traced)) return(eager)

  function(theta) {
    fn <- traced(theta)
    if (is.null(fn)) return(eager(theta))
    tt <- torch::torch_tensor(theta, dtype = torch::torch_float())
    as.numeric(torch::with_no_grad(fn(tt)))
  }
}

#' Replay the MDN's i.i.d. density as TorchScript instead of driving it from R
#'
#' Every operation in [mdn_iid_blocks()] crosses from R into libtorch, and at
#' MCMC batch sizes that crossing costs far more than the arithmetic behind it:
#' about 0.2 ms each, thirty of them per evaluation, against a few hundred
#' microseconds of actual work. [torch::jit_trace()] records the same
#' computation once and replays it inside libtorch, so an evaluation costs one
#' crossing rather than thirty.
#'
#' It is the same code either way. Tracing runs the eager path and records what
#' it did, so there is no second implementation of the density to keep in step
#' with this one or with the Stan generator.
#'
#' Three things make this a shortcut rather than the path. Recording a trace
#' costs several evaluations' worth of time, so nothing is recorded until the
#' evaluator has been called `warmup` times and it is clear this is a loop
#' rather than a one-off; a single [log_lik()] call should not pay for a
#' compiler. A trace fixes the shapes it was recorded at, so there is one per
#' number of parameter rows, checked against the eager result before anything
#' uses it. And it is only worth recording when the whole observation set fits
#' in one chunk, since otherwise the graph unrolls the chunk loop. Failing any
#' of these is not an error: the caller falls back to the eager path, which is
#' why `NULL` is a perfectly good answer here.
#'
#' Set `options(neuralsbi.jit = FALSE)` to skip tracing entirely.
#'
#' @param de,xt,max_batch As in [mdn_iid_blocks()], with the observations
#'   already a tensor.
#' @param eager The evaluator to check each trace against, and to fall back to.
#' @param warmup Calls to serve eagerly before recording anything.
#' @return `function(theta)` returning a traced function for that many rows, or
#'   `NULL` when the eager path should be used. `NULL` if tracing is switched
#'   off.
#' @keywords internal
mdn_trace_cache <- function(de, xt, max_batch, eager, warmup = 4L) {
  if (!isTRUE(getOption("neuralsbi.jit", TRUE))) return(NULL)
  n_obs <- xt$shape[1]
  budget <- de$n_components * de$dim_theta
  cache <- new.env(parent = emptyenv())
  calls <- 0L

  function(theta) {
    calls <<- calls + 1L
    n_theta <- nrow(theta)
    if (calls <= warmup) return(NULL)
    if (mdn_chunk_size(n_theta, max_batch, budget) < n_obs) return(NULL)
    key <- as.character(n_theta)
    hit <- cache[[key]]
    if (!is.null(hit)) return(hit$fn)

    tt <- torch::torch_tensor(theta, dtype = torch::torch_float())
    fn <- tryCatch(mdn_trace(de, xt, tt), error = function(e) NULL)
    if (!is.null(fn)) {
      # A traced graph has no data-dependent branches, so agreeing once at this
      # shape is the whole guarantee. Disagreeing means a shape was baked in
      # somewhere it should not have been, and the trace is discarded.
      ok <- tryCatch(
        isTRUE(all.equal(as.numeric(torch::with_no_grad(fn(tt))), eager(theta),
                         tolerance = 1e-5)),
        error = function(e) FALSE)
      if (!ok) fn <- NULL
    }
    cache[[key]] <- list(fn = fn)
    fn
  }
}

#' Record one trace of the summed density at the shape of `tt`
#'
#' Tracing refuses a captured tensor that carries a gradient, and the fitted
#' weights all do, so they are switched off around the recording and switched
#' back afterwards. Nothing else about the network changes.
#' @keywords internal
mdn_trace <- function(de, xt, tt) {
  pars <- de$net$parameters
  had_grad <- vapply(pars, function(w) w$requires_grad, logical(1))
  on.exit({
    for (i in seq_along(pars)) if (had_grad[i]) pars[[i]]$requires_grad_(TRUE)
  }, add = TRUE)
  for (w in pars) w$requires_grad_(FALSE)

  torch::jit_trace(function(theta) {
    mix <- mdn_mixture(de, theta)
    mdn_chunk_lp(mix, xt, de$dim_theta)$to(dtype = torch::torch_float64())$sum(dim = 2)
  }, tt)
}

#' Walk the theta blocks of the cross product, calling `collect(idx, lp)`
#'
#' `lp` is the flat vector [de_log_prob()][density_estimator] returns for the
#' block, running theta-major: the whole observation set for the first
#' parameter, then for the second, and so on.
#' @keywords internal
cross_iid <- function(de, x, theta, max_batch, collect) {
  n_theta <- nrow(theta)
  n_obs <- nrow(x)
  per_block <- max(1L, floor(max_batch / max(n_obs, 1L)))
  # The observation half of the cross product is the same in every block, so it
  # is built once at the widest one and truncated for a short final block.
  target <- x[rep(seq_len(n_obs), times = min(per_block, n_theta)), ,
              drop = FALSE]
  for (s in seq.int(1L, n_theta, by = per_block)) {
    idx <- seq.int(s, min(s + per_block - 1L, n_theta))
    rows <- seq_len(length(idx) * n_obs)
    collect(idx, de_log_prob(de, target[rows, , drop = FALSE],
                             theta[rep(idx, each = n_obs), , drop = FALSE]))
  }
  invisible(NULL)
}

#' Walk the observation chunks of an MDN's i.i.d. density, calling
#' `collect(idx, lp)` with the `(n_theta, n_chunk)` tensor of log-densities
#'
#' The MLP maps `theta` to the mixture parameters and never sees `x`, so the
#' forward pass and the Cholesky assembly run once for the whole observation
#' set and only the quadratic form is chunked.
#' @keywords internal
mdn_iid_blocks <- function(de, xt, theta, max_batch, collect) {
  n_theta <- nrow(theta)
  n_obs <- xt$shape[1]
  tt <- torch::torch_tensor(theta, dtype = torch::torch_float())
  torch::with_no_grad({
    mix <- mdn_mixture(de, tt)
    per_chunk <- mdn_chunk_size(n_theta, max_batch,
                                de$n_components * de$dim_theta)
    for (s in seq.int(1L, n_obs, by = per_chunk)) {
      idx <- seq.int(s, min(s + per_chunk - 1L, n_obs))
      # Slicing a tensor is not free; skip it when the chunk is everything.
      xs <- if (length(idx) == n_obs) xt else xt[idx, , drop = FALSE]
      collect(idx, mdn_chunk_lp(mix, xs, de$dim_theta))
    }
  })
  invisible(NULL)
}

#' The mixture an MDN puts over its target for each row of `tt`
#'
#' Split out from [mdn_iid_blocks()] so the eager driver and the traced one
#' share a single definition of the density. Everything here depends on the
#' conditioning variable alone, so it is computed once per call and reused
#' across every observation chunk.
#' @keywords internal
mdn_mixture <- function(de, tt) {
  params <- de$net(tt)
  L <- mdn_build_tril(de$net, params$tril_flat)                    # (T,K,p,p)
  diag_L <- torch::torch_diagonal(L, dim1 = 3, dim2 = 4)
  list(
    chol = L,
    log_w = torch::nnf_log_softmax(params$logits, dim = 2)$unsqueeze(3),
    logdet = 2 * torch::torch_log(diag_L)$sum(dim = 3)$unsqueeze(3),  # (T,K,1)
    means = params$means$unsqueeze(3)                                # (T,K,1,p)
  )
}

#' Log-density of a chunk of observations under each row's mixture
#'
#' @return A `(n_theta, n_chunk)` tensor.
#' @keywords internal
mdn_chunk_lp <- function(mix, xs, p) {
  # (T,K,n,p): every observation minus every component mean
  diff <- xs$unsqueeze(1)$unsqueeze(1) - mix$means
  z <- torch::linalg_solve_triangular(mix$chol, diff$transpose(3, 4),
                                      upper = FALSE)
  quad <- z$pow(2)$sum(dim = 3)                                    # (T,K,n)
  comp <- -0.5 * (p * log(2 * pi) + mix$logdet + quad)
  torch::torch_logsumexp(mix$log_w + comp, dim = 2)
}

#' How many observations to score in one call
#'
#' `max_batch` counts (theta, x) pairs, so that is what sizes the chunk. The
#' quadratic form materializes `K` components and `p` dimensions per pair on top
#' of that, so a second cap keeps the `(T, K, n, p)` intermediate to a few tens
#' of megabytes however many components the mixture has. Dividing by `K * p`
#' alone, as this used to, chunked a 5000-observation call ten ways and paid the
#' fixed cost of a torch call ten times for a 400 KB intermediate.
#' @keywords internal
mdn_chunk_size <- function(n_theta, max_batch, per_pair) {
  max(1L, min(floor(max_batch / max(n_theta, 1L)),
              floor(4e6 / max(n_theta * per_pair, 1L))))
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

  # A chain calls this thousands of times with the same observation, so
  # everything that depends on the observation alone is settled here: the
  # standardization, its Jacobian, and whatever the estimator wants to
  # precompute. Only the parameter side is left for the call itself.
  x_z <- apply_standardizer(fit$std_x, x_obs)
  offset <- nrow(x_z) * standardizer_log_jac(fit$std_x)
  ev <- de_iid_evaluator(fit$de, x_z, max_batch = max_batch)

  function(theta) {
    theta <- as_theta_matrix(theta, fit$dim_theta)
    lp <- as.numeric(prior$log_prob(theta))
    if (bounded) lp[!within_support(prior, theta)] <- -Inf
    ok <- is.finite(lp)
    if (!any(ok)) return(lp)
    # Only evaluate the network where the prior gives the point any mass; that
    # is a large saving once a chain starts probing a bounded edge.
    theta_z <- apply_standardizer(fit$std_theta, theta[ok, , drop = FALSE])
    lp[ok] <- lp[ok] + ev(theta_z) + offset
    lp
  }
}
