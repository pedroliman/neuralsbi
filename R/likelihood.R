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
  surrogate_score(fit, theta, x, sum_iid, max_batch)
}

#' How a surrogate fit is scored
#'
#' [nle()] and [nre()] fits are interchangeable everywhere downstream of the
#' estimator, and they differ in exactly three things: which function produces
#' the `n_theta x n_obs` matrix of scores, which one produces its row sums
#' without building it, and whether standardizing `x` needs a
#' change-of-variables term. A density reported in the simulator's units needs
#' one; a ratio does not, since the Jacobian cancels between its numerator and
#' denominator (see [nre()]).
#'
#' Keeping the three together in one table is what lets [surrogate_score()] and
#' [surrogate_potential()] be plain functions rather than a generic each. It
#' also means the fact that a ratio has no Jacobian is written down once.
#'
#' @param fit An `nsbi_nle` or `nsbi_nre` fit.
#' @return A list with `matrix_fn`, `evaluator` and `log_jac`.
#' @keywords internal
surrogate_ops <- function(fit) UseMethod("surrogate_ops")

#' @export
surrogate_ops.nsbi_nle <- function(fit) {
  # de_log_prob() works in standardized x space; the Jacobian puts it back in
  # the units the simulator returned.
  list(matrix_fn = de_log_lik_iid, evaluator = de_iid_evaluator,
       log_jac = TRUE)
}

#' The body [log_lik()] and [log_ratio()] share
#'
#' Both take a `theta` and an `x` whose rows are independent observations,
#' check them, standardize them, and either return the `n_theta x n_obs` matrix
#' or its row sums. Which functions do the scoring, and whether the Jacobian
#' comes with them, is [surrogate_ops()]'s business.
#'
#' @param fit An `nsbi_nle` or `nsbi_nre` fit.
#' @param theta,x,sum_iid,max_batch As in [log_lik()].
#' @keywords internal
surrogate_score <- function(fit, theta, x, sum_iid, max_batch) {
  check_fit_alive(fit)
  ops <- surrogate_ops(fit)
  matrix_fn <- ops$matrix_fn
  evaluator <- ops$evaluator
  log_jac <- ops$log_jac
  # check_matrix() only enforces type and shape, so an NA/NaN/Inf entry used to
  # standardize into another NA and come back as a silent NA log-lik or
  # log-ratio instead of an error (#202). Checked explicitly here, as
  # posterior.nsbi_npe() and mcmc_posterior() already do for the same fits.
  theta <- check_numeric(theta, "theta")
  check_finite(theta, "theta")
  x <- check_numeric(x, "x")
  check_finite(x, "x")

  # Both arguments are matrices with a required width, so a bare "expected 2
  # columns" would leave the user guessing which one is wrong.
  theta <- check_matrix(theta, fit$dim_theta, "theta",
                        "one parameter per column")
  x <- check_matrix(x, fit$dim_x, "x", "one row per independent observation")

  theta_z <- apply_standardizer(fit$std_theta, theta)
  x_z <- apply_standardizer(fit$std_x, x)

  # The Jacobian is constant, so the sum over observations picks it up once per
  # observation and the per-observation matrix picks it up in every cell.
  jac <- if (isTRUE(log_jac)) standardizer_log_jac(fit$std_x) else 0

  if (!isTRUE(sum_iid)) {
    return(matrix_fn(fit$de, x_z, theta_z, max_batch = max_batch) + jac)
  }
  ev <- evaluator(fit$de, x_z, max_batch = max_batch)
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
  iid_matrix(de, x, theta, max_batch, de_log_prob)
}

#' The `n_theta x n_obs` cross product, block by block
#'
#' The body of [de_log_lik_iid()]'s default method, with the per-pair scorer
#' left open so [nre()] can reuse it. `score(de, x_rows, theta_rows)` is
#' [de_log_prob()][density_estimator] for a density estimator and
#' [nre_score()] for a ratio estimator.
#' @inheritParams de_log_lik_iid
#' @param score The per-pair scorer, `function(de, x, theta)`.
#' @keywords internal
iid_matrix <- function(de, x, theta, max_batch, score) {
  n_obs <- nrow(x)
  out <- matrix(0, nrow = nrow(theta), ncol = n_obs)
  cross_iid(de, x, theta, max_batch, function(idx, lp) {
    out[idx, ] <<- matrix(lp, nrow = length(idx), ncol = n_obs, byrow = TRUE)
  }, score = score)
  out
}

#' @export
de_log_lik_iid.nsbi_de_lingauss <- function(de, x, theta, max_batch = 1e5) {
  # One conditional mean per parameter, then every observation scored against
  # it. No loop over observations at all.
  mu <- lingauss_mean(de, theta)
  lp <- vapply(seq_len(nrow(mu)),
               function(i) dmvnorm_chol(x, mu[i, ], de$chol),
               numeric(nrow(x)))
  # vapply returns a vector, not a one-row matrix, when there is a single
  # observation, so the shape is set explicitly rather than by transposing.
  matrix(lp, nrow = nrow(mu), ncol = nrow(x), byrow = TRUE)
}

#' @export
de_log_lik_iid.nsbi_de_mdn <- function(de, x, theta, max_batch = 1e5) {
  xt <- torch::torch_tensor(x, dtype = torch::torch_float(),
                            device = net_device(de$net))
  pieces <- list()
  mdn_iid_blocks(de, xt, theta, max_batch, function(idx, lp) {
    pieces[[length(pieces) + 1L]] <<-
      torch::as_array(lp$to(device = "cpu", dtype = torch::torch_float64()))
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
  iid_evaluator(de, x, max_batch, de_log_prob)
}

#' Summed log densities over a fixed observation set, block by block
#'
#' The body of [de_iid_evaluator()]'s default method, with the per-pair scorer
#' left open so [nre_iid_evaluator()] can reuse it.
#' @inheritParams de_iid_evaluator
#' @param score The per-pair scorer, `function(de, x, theta)`.
#' @keywords internal
iid_evaluator <- function(de, x, max_batch, score) {
  n_obs <- nrow(x)
  force(max_batch)
  force(score)
  function(theta) {
    out <- numeric(nrow(theta))
    cross_iid(de, x, theta, max_batch, function(idx, lp) {
      # `lp` runs theta-major, so a matrix of it has one column per theta.
      out[idx] <<- colSums(matrix(lp, nrow = n_obs))
    }, score = score)
    out
  }
}

#' @export
de_iid_evaluator.nsbi_de_lingauss <- function(de, x, max_batch = 1e5) {
  function(theta) {
    mu <- lingauss_mean(de, theta)
    vapply(seq_len(nrow(mu)),
           function(i) sum(dmvnorm_chol(x, mu[i, ], de$chol)),
           numeric(1))
  }
}

#' @export
de_iid_evaluator.nsbi_de_mdn <- function(de, x, max_batch = 1e5) {
  # The observations become a tensor once, not once per MCMC step.
  dev <- net_device(de$net)
  xt <- torch::torch_tensor(x, dtype = torch::torch_float(), device = dev)
  force(max_batch)
  eager <- function(theta) {
    total <- numeric(nrow(theta))
    mdn_iid_blocks(de, xt, theta, max_batch, function(idx, lp) {
      total <<- total +
        as.numeric(lp$to(device = "cpu",
                         dtype = torch::torch_float64())$sum(dim = 2))
    })
    total
  }
  traced <- mdn_trace_cache(de, xt, max_batch, eager)
  if (is.null(traced)) return(eager)

  function(theta) {
    fn <- traced(theta)
    if (is.null(fn)) return(eager(theta))
    tt <- torch::torch_tensor(theta, dtype = torch::torch_float(), device = dev)
    as.numeric(torch::with_no_grad(fn(tt))$to(device = "cpu"))
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
#' it did, so there is no separate implementation to keep in step with the
#' eager one -- that claim holds only for the traced path relative to
#' [mdn_log_prob_tensor()]. The MDN density still has three implementations in
#' total, one per runtime: [mdn_log_prob_tensor()] (the eager training path),
#' [mdn_mixture()]/[mdn_chunk_lp()] (the i.i.d. fast path used here), and
#' `stan_fn_mdn()` (the generated Stan code, `R/stan.R`). That is by design,
#' not drift: three runtimes need three implementations, and the tests pin
#' them to each other numerically.
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
#' @param warmup Calls to serve eagerly before recording anything. No caller
#'   overrides the default today, but it is a real tuning knob -- how many
#'   evaluations tracing costs before it pays for itself -- that a future
#'   caller would plausibly want to change, so it stays a parameter.
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

    tt <- torch::torch_tensor(theta, dtype = torch::torch_float(),
                              device = xt$device)
    fn <- tryCatch(mdn_trace(de, xt, tt), error = function(e) NULL)
    if (!is.null(fn)) {
      # A traced graph has no data-dependent branches, so agreeing once at this
      # shape is the whole guarantee. Disagreeing means a shape was baked in
      # somewhere it should not have been, and the trace is discarded.
      ok <- tryCatch(
        isTRUE(all.equal(
          as.numeric(torch::with_no_grad(fn(tt))$to(device = "cpu")),
          eager(theta), tolerance = 1e-5)),
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
#' `lp` is the flat vector `score()` returns for the block, running theta-major:
#' the whole observation set for the first parameter, then for the second, and
#' so on. `score` is [de_log_prob()][density_estimator] for a density estimator
#' and [nre_score()] for a ratio estimator; the blocking is the same either way.
#' @keywords internal
cross_iid <- function(de, x, theta, max_batch, collect, score = de_log_prob) {
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
    collect(idx, score(de, target[rows, , drop = FALSE],
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
  tt <- torch::torch_tensor(theta, dtype = torch::torch_float(),
                            device = xt$device)
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

#' Unnormalized log posterior of a surrogate fit
#'
#' \eqn{\log q_\phi(x \mid \theta) + \log p(\theta)} for an [nle()] fit and
#' \eqn{\log r_\phi(\theta, x) + \log p(\theta)} for an [nre()] fit, returning
#' `-Inf` outside the prior support. This is the potential the MCMC samplers
#' target, and the two fits differ only in what [surrogate_ops()] hands back.
#'
#' @param fit An `nsbi_nle` or `nsbi_nre` fit.
#' @param x_obs The observation to condition on. Rows are independent
#'   observations of the same parameter.
#' @param max_batch Largest number of `(theta, x)` pairs evaluated at once. No
#'   caller overrides it today, but it is a real tuning knob -- the batch size
#'   the MCMC evaluations are chunked into -- that a future caller would
#'   plausibly want to change, so it stays a parameter.
#' @return `function(theta)` giving one unnormalized log posterior density per
#'   row.
#' @keywords internal
surrogate_potential <- function(fit, x_obs, max_batch = 1e5) {
  ops <- surrogate_ops(fit)
  evaluator <- ops$evaluator
  log_jac <- ops$log_jac
  prior <- fit$prior
  # prior_custom() without a log_prob_fn returns NA rather than nothing, so the
  # only way to find out is to ask it. Better here than as a puzzling
  # initialization failure a few hundred lines later.
  probe <- tryCatch(prior$log_prob(sample_prior(prior, 2L)),
                    error = function(e) NA_real_)
  if (is.null(prior$log_prob) || all(is.na(probe))) {
    stop("Sampling this posterior needs a prior log-density, and this prior ",
         "does not have one.\nRebuild it with prior_custom(..., log_prob_fn = ).",
         call. = FALSE)
  }
  x_obs <- as_theta_matrix(x_obs, fit$dim_x)
  bounded <- !is.null(prior$lower) || !is.null(prior$upper)

  # A chain calls this thousands of times with the same observation, so
  # everything that depends on the observation alone is settled here: the
  # standardization, any Jacobian, and whatever the estimator wants to
  # precompute. Only the parameter side is left for the call itself.
  x_z <- apply_standardizer(fit$std_x, x_obs)
  offset <- if (isTRUE(log_jac)) {
    nrow(x_z) * standardizer_log_jac(fit$std_x)
  } else {
    0
  }
  ev <- evaluator(fit$de, x_z, max_batch = max_batch)

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
