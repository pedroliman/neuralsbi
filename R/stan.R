#' Export a learned likelihood to Stan
#'
#' `stan_code()` writes the likelihood learned by [nle()] out as Stan source: a
#' `functions` block that recomputes \eqn{\log q_\phi(x \mid \theta)} in Stan's
#' own language, with the trained weights travelling alongside as data. Stan
#' differentiates the generated code itself, so NUTS gets exact gradients and
#' nothing has to be linked against `torch` at run time.
#'
#' @section Why this is useful:
#'
#' The exported function is an ordinary Stan function of `theta`. That means
#' the surrogate likelihood stops being the whole model and becomes one term in
#' a model you write: put a hierarchical prior over the parameters, add
#' covariates, mix in a second data source with a likelihood you *do* know, and
#' let NUTS sample the lot. None of that is reachable from a posterior
#' estimator, which only ever knows the one conditional it was trained on.
#'
#' @section The generated interface:
#'
#' Two entry points, where `w` is the packed weight vector from [stan_data()]:
#'
#' ```stan
#' real nsbi_log_lik_lpdf(vector x, vector theta, vector w);      // one observation
#' real nsbi_log_lik_sum_lpdf(matrix x, vector theta, vector w);  // rows are i.i.d.
#' ```
#'
#' Both take `x` and `theta` in the original units the simulator and prior use;
#' the standardization the estimator trained under is folded into the generated
#' code. Prefer the `_sum` form for repeated observations: for the MDN and the
#' linear-Gaussian estimator the conditional distribution depends on `theta`
#' alone, so it is computed once and reused across every row.
#'
#' @section Why the code is generated rather than called:
#'
#' Python `sbi` hands its learned likelihood to Pyro or PyMC as a callable,
#' because everything lives in one process and the sampler can differentiate
#' the same graph the estimator was trained in. Stan cannot work that way: it
#' compiles to C++ and needs the gradient inside its own autodiff. The
#' alternative would be an external C++ header (`--allow-undefined` plus a
#' `USER_HEADER`) linking `torch` into every Stan compile, with the gradients
#' plumbed by hand. Generating source keeps the run-time dependency at zero and
#' makes the result something you can read, edit, and check -- which is what the
#' package's own tests do, by evaluating the emitted functions and comparing
#' them against [log_lik()].
#'
#' @section Supported estimators:
#'
#' `"mdn"`, `"maf"` and `"linear_gaussian"`. `"nsf"` is not exported -- its
#' rational-quadratic spline transform would be a large and fragile block of
#' generated Stan. Refit with `"maf"` if you need a flow.
#'
#' @param fit An `nsbi_nle` object from [nle()].
#' @param name Prefix for the generated functions.
#' @param model Generate a complete, runnable model (the default) or only the
#'   `functions` block, for `#include`-ing into a model of your own.
#' @param file Path to write to.
#' @param x_obs Observation to put in the data list. Rows are independent
#'   observations.
#'
#' @return `stan_code()` returns the Stan program as a single string;
#'   `write_stan_model()` returns `file` invisibly; `stan_data()` returns a
#'   named list ready for `cmdstanr` or `rstan`.
#'
#' @examples
#' prior <- prior_uniform(c(mu = -3), c(mu = 3))
#' fit <- nle(prior, function(mu) c(y = rnorm(1, mu, 0.5)),
#'            n_simulations = 500, density_estimator = "linear_gaussian")
#'
#' cat(substr(stan_code(fit), 1, 400))
#' str(stan_data(fit, matrix(rnorm(10), ncol = 1)), max.level = 1)
#' @name stan_export
NULL

#' @rdname stan_export
#' @export
stan_code <- function(fit, name = "nsbi_log_lik", model = TRUE) {
  stopifnot(inherits(fit, "nsbi_nle"))
  check_fit_alive(fit)
  if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", name)) {
    stop("`name` must be a valid Stan identifier.", call. = FALSE)
  }
  packed <- stan_pack(fit)
  fns <- stan_functions_block(fit, name, packed)
  if (!isTRUE(model)) return(fns)
  paste0(fns, "\n", stan_model_blocks(fit, name, packed))
}

#' @rdname stan_export
#' @export
write_stan_model <- function(fit, file, name = "nsbi_log_lik", model = TRUE) {
  # Before the code is generated: transpiling the network only to fail on the
  # destination is work thrown away, and writeLines() takes a connection as
  # well as a path, so a wrong value gets some distance before it complains.
  check_path(file, "file")
  writeLines(stan_code(fit, name = name, model = model), file)
  invisible(file)
}

#' @rdname stan_export
#' @export
stan_data <- function(fit, x_obs = NULL) {
  stopifnot(inherits(fit, "nsbi_nle"))
  # Both halves of the export read the weights, so both need a live network.
  # Without this, a fit restored by readRDS() gets the "save with save_npe()"
  # message from stan_code() and a dangling-pointer error from net_param()
  # here, for the same fit and the same cause.
  check_fit_alive(fit)
  packed <- stan_pack(fit)
  out <- list(nsbi_nw = length(packed$w), nsbi_w = packed$w)
  if (!is.null(x_obs)) {
    x_obs <- as_theta_matrix(x_obs, fit$dim_x)
    out$N <- nrow(x_obs)
    out$x <- x_obs
    dimnames(out$x) <- NULL
  }
  prior <- fit$prior
  if (identical(prior$type, "uniform")) {
    out$nsbi_low <- as.numeric(prior$lower)
    out$nsbi_high <- as.numeric(prior$upper)
  } else if (identical(prior$type, "normal")) {
    out$nsbi_prior_mean <- as.numeric(prior$params$mean)
    out$nsbi_prior_sd <- as.numeric(prior$params$sd)
  }
  out
}

# ---- weight packing -------------------------------------------------------

#' Flatten the trained weights into one vector, recording where each block sits
#'
#' Stan has no closures and no struct types, so everything the generated
#' function needs arrives as a single `vector w` and is sliced back out inside.
#' Matrices are flattened column-major, which is exactly what Stan's
#' `to_matrix()` expects, so the two sides agree without any transposition.
#' @keywords internal
stan_pack <- function(fit) {
  p <- new.env(parent = emptyenv())
  p$parts <- list()
  p$offset <- 1L
  p$blocks <- list()

  add <- function(key, value, nrow = NULL, ncol = NULL) {
    v <- as.numeric(value)
    p$parts[[length(p$parts) + 1L]] <- v
    p$blocks[[key]] <- list(offset = p$offset, n = length(v),
                            nrow = nrow, ncol = ncol)
    p$offset <- p$offset + length(v)
    invisible(NULL)
  }

  de <- fit$de
  kind <- class(de)[1L]
  switch(
    kind,
    nsbi_de_lingauss = stan_pack_lingauss(de, add),
    nsbi_de_mdn = stan_pack_mdn(de, add),
    nsbi_de_maf = stan_pack_maf(de, add),
    stop(sprintf(
      paste0("Cannot export a '%s' estimator to Stan.\n",
             "Supported: \"mdn\", \"maf\", \"linear_gaussian\". The NSF's ",
             "spline transform is not generated; refit with \"maf\"."),
      fit$density_estimator), call. = FALSE)
  )

  list(w = unlist(p$parts, use.names = FALSE), blocks = p$blocks, kind = kind)
}

#' A torch tensor as a plain R matrix, in double precision
#' @keywords internal
torch_mat <- function(t) {
  m <- torch::as_array(t$to(device = "cpu", dtype = torch::torch_float64()))
  if (is.null(dim(m))) matrix(m, ncol = 1L) else m
}

#' Named parameter lookup that fails loudly rather than silently returning NULL
#' @keywords internal
net_param <- function(net, name) {
  params <- net$parameters
  if (!name %in% names(params)) {
    stop(sprintf("Expected parameter '%s' on the network; found: %s.",
                 name, paste(names(params), collapse = ", ")), call. = FALSE)
  }
  torch_mat(params[[name]])
}

#' @keywords internal
stan_pack_lingauss <- function(de, add) {
  # B is (1 + dim_cond) x dim_target; chol() is the upper factor, Stan wants
  # the lower one.
  add("B", de$B, nrow = nrow(de$B), ncol = ncol(de$B))
  add("L", t(de$chol), nrow = nrow(de$chol), ncol = ncol(de$chol))
}

#' @keywords internal
stan_pack_mdn <- function(de, add) {
  require_torch()
  net <- de$net
  for (l in seq_along(de$hidden)) {
    W <- net_param(net, sprintf("trunk.%d.weight", 2L * (l - 1L)))
    b <- net_param(net, sprintf("trunk.%d.bias", 2L * (l - 1L)))
    add(sprintf("W%d", l), W, nrow = nrow(W), ncol = ncol(W))
    add(sprintf("b%d", l), b, nrow = length(b), ncol = 1L)
  }
  for (head in c("logits", "means", "tril")) {
    W <- net_param(net, sprintf("head_%s.weight", head))
    b <- net_param(net, sprintf("head_%s.bias", head))
    add(paste0("W_", head), W, nrow = nrow(W), ncol = ncol(W))
    add(paste0("b_", head), b, nrow = length(b), ncol = 1L)
  }
}

#' @keywords internal
stan_pack_maf <- function(de, add) {
  require_torch()
  net <- de$net
  # Masks are recomputed rather than read off the module: they are a pure
  # function of the shapes, and folding them into the exported weights means
  # the generated Stan sees ordinary dense matrices.
  masks <- made_masks(de$dim_theta, de$dim_x, de$hidden)
  for (k in seq_len(de$n_transforms)) {
    pre <- sprintf("mades.%d", k - 1L)   # nn_module_list names are 0-indexed
    for (l in seq_along(de$hidden)) {
      W <- net_param(net, sprintf("%s.trunk.%d.linear.weight", pre, 2L * (l - 1L)))
      b <- net_param(net, sprintf("%s.trunk.%d.linear.bias", pre, 2L * (l - 1L)))
      add(sprintf("W%d_%d", k, l), W * masks$hidden[[l]],
          nrow = nrow(W), ncol = ncol(W))
      add(sprintf("b%d_%d", k, l), b, nrow = length(b), ncol = 1L)
    }
    for (head in c("mu", "alpha")) {
      W <- net_param(net, sprintf("%s.out_%s.linear.weight", pre, head))
      b <- net_param(net, sprintf("%s.out_%s.linear.bias", pre, head))
      add(sprintf("W%s%d", head, k), W * masks$out, nrow = nrow(W), ncol = ncol(W))
      add(sprintf("b%s%d", head, k), b, nrow = length(b), ncol = 1L)
    }
  }
}

# ---- code generation ------------------------------------------------------

#' @keywords internal
stan_num <- function(x) {
  vapply(as.numeric(x), function(v) sprintf("%.17g", v), character(1))
}

#' A Stan vector literal, e.g. `[1, 2, 3]'`
#' @keywords internal
stan_vec <- function(x) {
  paste0("[", paste(stan_num(x), collapse = ", "), "]'")
}

#' Stan expression slicing a matrix block back out of `w`
#' @keywords internal
stan_mat_of <- function(blocks, key, wname = "w") {
  b <- blocks[[key]]
  sprintf("to_matrix(segment(%s, %d, %d), %d, %d)",
          wname, b$offset, b$n, b$nrow, b$ncol)
}

#' Stan expression slicing a vector block back out of `w`
#' @keywords internal
stan_vec_of <- function(blocks, key, wname = "w") {
  b <- blocks[[key]]
  sprintf("segment(%s, %d, %d)", wname, b$offset, b$n)
}

#' @keywords internal
stan_functions_block <- function(fit, name, packed) {
  body <- switch(
    packed$kind,
    nsbi_de_lingauss = stan_fn_lingauss(fit, name, packed),
    nsbi_de_mdn = stan_fn_mdn(fit, name, packed),
    nsbi_de_maf = stan_fn_maf(fit, name, packed)
  )
  paste0(
    "// Generated by neuralsbi::stan_code(). Do not edit by hand.\n",
    sprintf("// Surrogate likelihood q(x | theta) from an nle() fit (%s, %d simulations).\n",
            fit$density_estimator, fit$n_simulations),
    "functions {\n",
    if (packed$kind == "nsbi_de_lingauss") "" else stan_relu_helper(name),
    body,
    "}\n"
  )
}

#' `+ c` or `- |c|`, so the generated arithmetic reads like arithmetic
#' @keywords internal
stan_addend <- function(x) {
  x <- as.numeric(x)
  sprintf("%s %s", if (x < 0) "-" else "+", stan_num(abs(x)))
}

#' @keywords internal
stan_relu_helper <- function(name) {
  paste0(
    sprintf("  vector %s_relu(vector v) {\n", name),
    "    vector[num_elements(v)] out;\n",
    "    for (i in 1:num_elements(v)) out[i] = fmax(v[i], 0.0);\n",
    "    return out;\n",
    "  }\n\n"
  )
}

#' Standardization preamble shared by every generated entry point
#' @keywords internal
stan_standardize_lines <- function(fit, target = "x", cond = "theta") {
  paste0(
    sprintf("    vector[%d] xs = (%s - %s) ./ %s;\n",
            fit$dim_x, target, stan_vec(fit$std_x$center), stan_vec(fit$std_x$scale)),
    sprintf("    vector[%d] ts = (%s - %s) ./ %s;\n",
            fit$dim_theta, cond, stan_vec(fit$std_theta$center),
            stan_vec(fit$std_theta$scale))
  )
}

#' The constant putting a standardized density back into the data's own units
#' @keywords internal
stan_jacobian <- function(fit) stan_num(standardizer_log_jac(fit$std_x))

#' The i.i.d.-sum entry point shared by every generated `_sum_lpdf`
#'
#' Standardizes `theta`, declares `x`'s center/scale, and accumulates `body`
#' (an expression for one observation's log density, in terms of `x[n]`) over
#' `rows(x)` before applying the jacobian once at the end. `precompute` is the
#' one place estimators differ: `linear_gaussian` and the MDN can build their
#' conditional distribution once, outside the loop, because it depends on
#' `theta` alone; MAF has nothing to hoist, so its caller leaves this at the
#' default. Mirrors [de_log_lik_iid()] (R/likelihood.R), which sums the same
#' per-observation log density on the R side; the two have to agree.
#' @keywords internal
stan_sum_lines <- function(fit, P, body, precompute = "") {
  paste0(
    sprintf("    vector[%d] ts = (theta - %s) ./ %s;\n",
            fit$dim_theta, stan_vec(fit$std_theta$center),
            stan_vec(fit$std_theta$scale)),
    precompute,
    sprintf("    vector[%d] xc = %s;\n", P, stan_vec(fit$std_x$center)),
    sprintf("    vector[%d] xsc = %s;\n", P, stan_vec(fit$std_x$scale)),
    "    real total = 0;\n",
    "    for (n in 1:rows(x)) {\n",
    sprintf("      total += %s;\n", body),
    "    }\n",
    sprintf("    return total + rows(x) * (%s);\n", stan_jacobian(fit))
  )
}

# ---- linear-Gaussian ------------------------------------------------------

#' @keywords internal
stan_fn_lingauss <- function(fit, name, packed) {
  b <- packed$blocks
  P <- fit$dim_x
  paste0(
    sprintf("  // conditional mean, given the standardized parameter\n"),
    sprintf("  vector %s_mean(vector ts, vector w) {\n", name),
    sprintf("    return (%s)' * append_row(1.0, ts);\n", stan_mat_of(b, "B")),
    "  }\n\n",
    sprintf("  real %s_lpdf(vector x, vector theta, vector w) {\n", name),
    stan_standardize_lines(fit),
    sprintf("    return multi_normal_cholesky_lpdf(xs | %s_mean(ts, w), %s) %s;\n",
            name, stan_mat_of(b, "L"), stan_addend(standardizer_log_jac(fit$std_x))),
    "  }\n\n",
    sprintf("  real %s_sum_lpdf(matrix x, vector theta, vector w) {\n", name),
    stan_sum_lines(
      fit, P,
      body = "multi_normal_cholesky_lpdf((x[n]' - xc) ./ xsc | mu, L)",
      # The conditional distribution depends on theta alone, so it is built
      # once and reused for every observation.
      precompute = paste0(
        sprintf("    vector[%d] mu = %s_mean(ts, w);\n", P, name),
        sprintf("    matrix[%d, %d] L = %s;\n", P, P, stan_mat_of(b, "L"))
      )
    ),
    "  }\n"
  )
}

# ---- MDN ------------------------------------------------------------------

#' @keywords internal
stan_fn_mdn <- function(fit, name, packed) {
  b <- packed$blocks
  de <- fit$de
  P <- de$dim_theta          # the estimator's target: the data dimension
  K <- de$n_components
  tri <- tril_indices(P)
  Tn <- length(tri$row)
  head_len <- K + K * P + K * Tn

  trunk <- ""
  prev <- "ts"
  for (l in seq_along(de$hidden)) {
    trunk <- paste0(trunk, sprintf(
      "    vector[%d] h%d = %s_relu(%s * %s + %s);\n",
      de$hidden[l], l, name, stan_mat_of(b, sprintf("W%d", l)), prev,
      stan_vec_of(b, sprintf("b%d", l))))
    prev <- sprintf("h%d", l)
  }

  # Rebuilding one component's Cholesky factor from the flat lower-triangular
  # head, matching mdn_build_tril(): softplus on the diagonal, raw elsewhere.
  tril_lines <- vapply(seq_len(Tn), function(m) {
    src <- sprintf("tflat[(k - 1) * %d + %d]", Tn, m)
    val <- if (tri$is_diag[m]) sprintf("log1p_exp(%s) + 1e-6", src) else src
    sprintf("        L[%d, %d] = %s;\n", tri$row[m], tri$col[m], val)
  }, character(1))

  paste0(
    "  // MLP head: mixture logits, means and Cholesky factors, all functions\n",
    "  // of theta alone. Returned packed so the i.i.d. sum can hoist it out\n",
    "  // of the observation loop.\n",
    sprintf("  vector %s_head(vector ts, vector w) {\n", name),
    trunk,
    sprintf("    return append_row(append_row(%s * %s + %s, %s * %s + %s), %s * %s + %s);\n",
            stan_mat_of(b, "W_logits"), prev, stan_vec_of(b, "b_logits"),
            stan_mat_of(b, "W_means"), prev, stan_vec_of(b, "b_means"),
            stan_mat_of(b, "W_tril"), prev, stan_vec_of(b, "b_tril")),
    "  }\n\n",
    sprintf("  real %s_from_head(vector xs, vector head) {\n", name),
    sprintf("    vector[%d] logits = head[1:%d];\n", K, K),
    sprintf("    vector[%d] mflat = head[%d:%d];\n", K * P, K + 1L, K + K * P),
    sprintf("    vector[%d] tflat = head[%d:%d];\n", K * Tn, K + K * P + 1L, head_len),
    sprintf("    vector[%d] lp;\n", K),
    sprintf("    for (k in 1:%d) {\n", K),
    sprintf("      matrix[%d, %d] L = rep_matrix(0.0, %d, %d);\n", P, P, P, P),
    sprintf("      vector[%d] mu = mflat[((k - 1) * %d + 1):(k * %d)];\n", P, P, P),
    "      {\n", paste(tril_lines, collapse = ""), "      }\n",
    "      lp[k] = logits[k] + multi_normal_cholesky_lpdf(xs | mu, L);\n",
    "    }\n",
    "    return log_sum_exp(lp) - log_sum_exp(logits);\n",
    "  }\n\n",
    sprintf("  real %s_lpdf(vector x, vector theta, vector w) {\n", name),
    stan_standardize_lines(fit),
    sprintf("    return %s_from_head(xs, %s_head(ts, w)) %s;\n",
            name, name, stan_addend(standardizer_log_jac(fit$std_x))),
    "  }\n\n",
    sprintf("  real %s_sum_lpdf(matrix x, vector theta, vector w) {\n", name),
    stan_sum_lines(
      fit, P,
      body = sprintf("%s_from_head((x[n]' - xc) ./ xsc, head)", name),
      precompute = sprintf("    vector[%d] head = %s_head(ts, w);\n", head_len, name)
    ),
    "  }\n"
  )
}

# ---- MAF ------------------------------------------------------------------

#' @keywords internal
stan_fn_maf <- function(fit, name, packed) {
  b <- packed$blocks
  de <- fit$de
  P <- de$dim_theta          # the estimator's target: the data dimension
  rev_idx <- rev(seq_len(P))

  steps <- ""
  for (k in seq_len(de$n_transforms)) {
    if (k > 1L && P > 1L) {
      steps <- paste0(steps, sprintf(
        "    z = z[{%s}];\n", paste(rev_idx, collapse = ", ")))
    }
    prev <- sprintf("append_row(z, ts)")
    for (l in seq_along(de$hidden)) {
      steps <- paste0(steps, sprintf(
        "    vector[%d] h%d_%d = %s_relu(%s * %s + %s);\n",
        de$hidden[l], k, l, name, stan_mat_of(b, sprintf("W%d_%d", k, l)),
        prev, stan_vec_of(b, sprintf("b%d_%d", k, l))))
      prev <- sprintf("h%d_%d", k, l)
    }
    steps <- paste0(steps, sprintf(
      "    vector[%d] mu%d = %s * %s + %s;\n",
      P, k, stan_mat_of(b, sprintf("Wmu%d", k)), prev,
      stan_vec_of(b, sprintf("bmu%d", k))))
    # The clamp mirrors made_module()'s torch_clamp on alpha; without it a
    # saturating transform would disagree with the R fit in the tails.
    steps <- paste0(steps, sprintf(
      "    vector[%d] al%d;\n", P, k), sprintf(
      "    for (i in 1:%d) al%d[i] = fmin(fmax((%s * %s + %s)[i], -8.0), 8.0);\n",
      P, k, stan_mat_of(b, sprintf("Walpha%d", k)), prev,
      stan_vec_of(b, sprintf("balpha%d", k))))
    steps <- paste0(steps, sprintf(
      "    z = (z - mu%d) .* exp(-al%d);\n    logdet -= sum(al%d);\n", k, k, k))
  }

  paste0(
    "  // One forward pass through the stacked MADE transforms. Unlike the MDN,\n",
    "  // every transform depends on the data as well as theta, so an i.i.d.\n",
    "  // sum has to repeat the whole pass per observation.\n",
    sprintf("  real %s_std(vector xs, vector ts, vector w) {\n", name),
    sprintf("    vector[%d] z = xs;\n", P),
    "    real logdet = 0;\n",
    steps,
    sprintf("    return -0.5 * (dot_self(z) + %d * log(2 * pi())) + logdet;\n", P),
    "  }\n\n",
    sprintf("  real %s_lpdf(vector x, vector theta, vector w) {\n", name),
    stan_standardize_lines(fit),
    sprintf("    return %s_std(xs, ts, w) %s;\n", name,
            stan_addend(standardizer_log_jac(fit$std_x))),
    "  }\n\n",
    sprintf("  real %s_sum_lpdf(matrix x, vector theta, vector w) {\n", name),
    stan_sum_lines(fit, P, body = sprintf("%s_std((x[n]' - xc) ./ xsc, ts, w)", name)),
    "  }\n"
  )
}

# ---- the surrounding model ------------------------------------------------

#' Data, parameters and model blocks around the generated likelihood
#'
#' Only `prior_uniform()` and `prior_normal()` can be written out; anything
#' built with `prior_custom()` is arbitrary R code with no Stan counterpart, so
#' the model block is the user's to write.
#' @keywords internal
stan_model_blocks <- function(fit, name, packed) {
  prior <- fit$prior
  Q <- fit$dim_theta
  P <- fit$dim_x

  decl <- switch(
    prior$type %||% "custom",
    uniform = sprintf("  vector<lower=nsbi_low, upper=nsbi_high>[%d] theta;\n", Q),
    normal = sprintf("  vector[%d] theta;\n", Q),
    stop("Only prior_uniform() and prior_normal() can be written out as Stan.\n",
         "For any other prior, take stan_code(fit, model = FALSE) and write ",
         "the model block yourself.", call. = FALSE)
  )
  prior_data <- switch(
    prior$type,
    uniform = sprintf("  vector[%d] nsbi_low;\n  vector[%d] nsbi_high;\n", Q, Q),
    normal = sprintf("  vector[%d] nsbi_prior_mean;\n  vector[%d] nsbi_prior_sd;\n", Q, Q)
  )
  prior_stmt <- switch(
    prior$type,
    # A uniform prior is implicit in the declared bounds.
    uniform = "",
    normal = "  theta ~ normal(nsbi_prior_mean, nsbi_prior_sd);\n"
  )

  paste0(
    "data {\n",
    "  int<lower=1> N;                 // number of independent observations\n",
    sprintf("  matrix[N, %d] x;                 // the observations\n", P),
    "  int<lower=1> nsbi_nw;\n",
    "  vector[nsbi_nw] nsbi_w;         // trained weights, from stan_data()\n",
    prior_data,
    "}\n\n",
    "parameters {\n", decl, "}\n\n",
    "model {\n",
    prior_stmt,
    sprintf("  x ~ %s_sum(theta, nsbi_w);\n", name),
    "}\n"
  )
}

# ---- running it -----------------------------------------------------------

#' Sample an NLE posterior with Stan
#'
#' Writes the model, compiles it, and runs NUTS. Prefers \pkg{cmdstanr} and
#' falls back to \pkg{rstan}.
#' @keywords internal
stan_sample_nle <- function(fit, x_obs, ctl, n, verbose = FALSE) {
  dots <- ctl$dots
  code <- stan_code(fit)
  data <- stan_data(fit, x_obs)

  # NUTS draws are close to independent already, so `thin` -- which exists to
  # decorrelate the slice sampler's output -- does not apply here.
  iter_sampling <- dots$iter_sampling %||% ceiling(n / ctl$n_chains)
  iter_warmup <- dots$iter_warmup %||% max(ctl$warmup, 200L)
  refresh <- dots$refresh %||% (if (isTRUE(verbose)) 100L else 0L)

  draws <- if (requireNamespace("cmdstanr", quietly = TRUE)) {
    stan_run_cmdstanr(code, data, ctl, iter_warmup, iter_sampling, refresh)
  } else if (requireNamespace("rstan", quietly = TRUE)) {
    stan_run_rstan(code, data, ctl, iter_warmup, iter_sampling, refresh)
  } else {
    stop("sampler = \"stan\" needs cmdstanr or rstan installed.\n",
         "Install one, or use the built-in sampler with sampler = \"slice\".",
         call. = FALSE)
  }

  # draws arrives as iterations x chains x dim
  list(draws = matrix(aperm(draws, c(2, 1, 3)), ncol = dim(draws)[3]),
       diagnostics = mcmc_diagnostics(draws))
}

#' @keywords internal
stan_run_cmdstanr <- function(code, data, ctl, iter_warmup, iter_sampling,
                              refresh) {
  f <- tempfile(fileext = ".stan")
  on.exit(unlink(f), add = TRUE)
  writeLines(code, f)
  mod <- cmdstanr::cmdstan_model(f)
  fit <- mod$sample(data = data, chains = ctl$n_chains,
                    parallel_chains = min(ctl$n_chains, stan_cores()),
                    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                    seed = ctl$seed %||% sample.int(.Machine$integer.max, 1L),
                    refresh = refresh, show_messages = refresh > 0L)
  d <- fit$draws(variables = "theta", format = "draws_array")
  array(as.numeric(d), dim = dim(d))
}

#' @keywords internal
stan_run_rstan <- function(code, data, ctl, iter_warmup, iter_sampling,
                           refresh) {
  mod <- rstan::stan_model(model_code = code)
  fit <- rstan::sampling(mod, data = data, chains = ctl$n_chains,
                         cores = min(ctl$n_chains, stan_cores()),
                         warmup = iter_warmup,
                         iter = iter_warmup + iter_sampling,
                         seed = ctl$seed %||% sample.int(.Machine$integer.max, 1L),
                         refresh = refresh)
  # permuted = FALSE keeps the chains apart, which is what the diagnostics
  # want, and puts them in the iterations x chains x dim order cmdstanr's
  # draws_array already uses.
  d <- rstan::extract(fit, pars = "theta", permuted = FALSE)
  array(as.numeric(d), dim = dim(d))
}

#' @keywords internal
stan_cores <- function() {
  # detectCores() returns NA, not NULL, when it cannot tell.
  n <- parallel::detectCores(logical = FALSE)
  if (!is.finite(n)) n <- 1L
  max(1L, min(as.integer(n), 4L))
}
