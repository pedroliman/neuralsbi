#' Neural Spline Flow (NSF) conditional density estimator
#'
#' An autoregressive flow whose per-dimension transform is a monotonic
#' rational-quadratic spline (Durkan et al., 2019) instead of MAF's affine
#' shift-and-scale. Splines are far more expressive per layer, which helps on
#' sharply non-Gaussian posteriors (SLCP, two moons). We reuse the MADE
#' masking machinery from `R/flows.R`; each MADE outputs `3K - 1` spline
#' parameters per dimension (K bin widths, K bin heights, K - 1 interior
#' derivatives). The spline acts on `[-B, B]` and is the identity outside
#' (linear tails), so the standard-normal base distribution is unaffected in
#' the tails. Note: NSF implementations elsewhere often use coupling layers;
#' ours is autoregressive — same density family, different conditioning
#' structure.
#'
#' Select with `npe(..., density_estimator = "nsf")`.
#'
#' @keywords internal
#' @name nsf
NULL

#' Monotonic rational-quadratic spline, batched.
#'
#' @param inputs `(N,)` tensor of values to transform.
#' @param w_un,h_un,d_un Unnormalized widths `(N, K)`, heights `(N, K)`, and
#'   interior derivatives `(N, K - 1)`.
#' @param inverse Apply the inverse transform.
#' @param tail_bound Spline acts on `[-tail_bound, tail_bound]`; identity
#'   outside.
#' @return `list(outputs, logdet)`, both `(N,)` tensors.
#' @keywords internal
rq_spline <- function(inputs, w_un, h_un, d_un, inverse = FALSE,
                      tail_bound = 3, min_bin = 1e-3, min_deriv = 1e-3) {
  K <- w_un$shape[2]
  B <- tail_bound
  inside <- (inputs >= -B) & (inputs <= B)
  outputs <- inputs$clone()
  logdet <- torch::torch_zeros_like(inputs)
  if (as.numeric(inside$sum()) == 0) {
    return(list(outputs = outputs, logdet = logdet))
  }
  idx <- torch::torch_nonzero(inside)$squeeze(2)
  xin <- inputs[idx]
  w_un <- w_un[idx, , drop = FALSE]
  h_un <- h_un[idx, , drop = FALSE]
  d_un <- d_un[idx, , drop = FALSE]

  widths <- torch::nnf_softmax(w_un, dim = 2)
  widths <- min_bin + (1 - min_bin * K) * widths
  heights <- torch::nnf_softmax(h_un, dim = 2)
  heights <- min_bin + (1 - min_bin * K) * heights
  cumw <- torch::torch_cumsum(widths, dim = 2) * (2 * B) - B
  cumh <- torch::torch_cumsum(heights, dim = 2) * (2 * B) - B
  pad <- torch::torch_full(c(cumw$shape[1], 1), -B)
  cumw <- torch::torch_cat(list(pad, cumw), dim = 2)   # (N, K+1)
  cumh <- torch::torch_cat(list(pad, cumh), dim = 2)
  # exact endpoints (cumsum can drift slightly)
  cumw[, K + 1] <- B
  cumh[, K + 1] <- B
  # derivatives: softplus, boundary derivatives fixed at 1 for linear tails
  derivs <- min_deriv + torch::nnf_softplus(d_un + 0.5413)  # ~1 at init 0
  ones <- torch::torch_ones(c(derivs$shape[1], 1))
  derivs <- torch::torch_cat(list(ones, derivs, ones), dim = 2)  # (N, K+1)

  bin_grid <- if (inverse) cumh else cumw
  # searchsorted(right = TRUE) counts elements <= input, which for a value in
  # bin k of the (K+1)-point grid is exactly k; clamp handles the endpoints.
  k <- torch::torch_searchsorted(bin_grid, xin$unsqueeze(2), right = TRUE)
  k <- k$clamp(min = 1L, max = K)                                # (N, 1)

  g <- function(t, i) t$gather(2, i)$squeeze(2)
  xk <- g(cumw, k);  xk1 <- g(cumw, k + 1L)
  yk <- g(cumh, k);  yk1 <- g(cumh, k + 1L)
  dk <- g(derivs, k); dk1 <- g(derivs, k + 1L)
  w <- xk1 - xk
  h <- yk1 - yk
  s <- h / w

  if (!inverse) {
    xi <- (xin - xk) / w
    xi1m <- 1 - xi
    num <- h * (s * xi$pow(2) + dk * xi * xi1m)
    den <- s + (dk1 + dk - 2 * s) * xi * xi1m
    yout <- yk + num / den
    dnum <- s$pow(2) * (dk1 * xi$pow(2) + 2 * s * xi * xi1m + dk * xi1m$pow(2))
    ld <- torch::torch_log(dnum) - 2 * torch::torch_log(den)
  } else {
    term <- xin - yk
    r <- dk1 + dk - 2 * s
    a <- h * (s - dk) + term * r
    b <- h * dk - term * r
    cc <- -s * term
    disc <- b$pow(2) - 4 * a * cc
    xi <- 2 * cc / (-b - torch::torch_sqrt(disc$clamp(min = 0)))
    xi <- xi$clamp(min = 0, max = 1)
    xi1m <- 1 - xi
    yout <- xk + xi * w
    den <- s + r * xi * xi1m
    dnum <- s$pow(2) * (dk1 * xi$pow(2) + 2 * s * xi * xi1m + dk * xi1m$pow(2))
    ld <- -(torch::torch_log(dnum) - 2 * torch::torch_log(den))
  }
  outputs <- torch::torch_index_put(outputs, list(idx), yout)
  logdet <- torch::torch_index_put(logdet, list(idx), ld)
  list(outputs = outputs, logdet = logdet)
}

#' MADE block emitting 3K - 1 spline parameters per dimension
#' @keywords internal
nsf_made_module <- function(dim_x, dim_theta, hidden, n_bins) {
  torch::nn_module(
    classname = "nsbi_nsf_made",
    initialize = function() {
      m <- made_masks(dim_theta, dim_x, hidden)
      sizes <- c(dim_theta + dim_x, hidden)
      layers <- list()
      for (l in seq_along(hidden)) {
        layers[[length(layers) + 1L]] <-
          masked_linear(sizes[l], sizes[l + 1L], m$hidden[[l]])
        layers[[length(layers) + 1L]] <- torch::nn_relu()
      }
      self$trunk <- do.call(torch::nn_sequential, layers)
      self$n_params <- 3L * n_bins - 1L
      self$n_bins <- n_bins
      self$dim_theta <- dim_theta
      # output (i, d, j) reads flat index (d-1)*n_params + j -> repeat rows
      out_mask <- m$out[rep(seq_len(dim_theta), each = self$n_params), ,
                        drop = FALSE]
      self$out <- masked_linear(sizes[length(sizes)],
                                dim_theta * self$n_params, out_mask)
      # near-zero init => near-identity spline at start of training
      torch::nn_init_normal_(self$out$linear$weight, std = 1e-3)
      torch::nn_init_zeros_(self$out$linear$bias)
      self$out$linear$weight$data()$mul_(self$out$mask)
    },
    forward = function(theta, x) {
      h <- self$trunk(torch::torch_cat(list(theta, x), dim = 2))
      params <- self$out(h)$view(c(-1, self$dim_theta, self$n_params))
      K <- self$n_bins
      list(w = params[, , 1:K],
           h = params[, , (K + 1):(2 * K)],
           d = params[, , (2 * K + 1):(3 * K - 1)])
    }
  )
}

#' The full NSF: stacked spline-autoregressive transforms with order reversal
#' @keywords internal
nsf_module <- function(dim_x, dim_theta, n_transforms, hidden, n_bins,
                       tail_bound, embedding = NULL) {
  torch::nn_module(
    classname = "nsbi_nsf_net",
    initialize = function() {
      self$dim_theta <- dim_theta
      self$n_transforms <- n_transforms
      self$n_bins <- n_bins
      self$tail_bound <- tail_bound
      self$has_embedding <- !is.null(embedding)
      if (self$has_embedding) {
        self$embedding <- build_embedding_module(embedding, dim_x)
      }
      cond_dim <- embedding_output_dim(embedding, dim_x)
      self$mades <- torch::nn_module_list(
        lapply(seq_len(n_transforms),
               function(k) nsf_made_module(cond_dim, dim_theta, hidden, n_bins)())
      )
    },
    forward = function(theta, x) {
      stop("call nsf_log_prob_tensor() / de_sample() instead")
    }
  )
}

#' Apply one spline transform elementwise over all theta dimensions.
#'
#' Spline parameters are computed from `z` (and `x`); the transform itself is
#' applied to `values` (defaults to `z`). The separation matters when
#' inverting: parameters must come from the partially reconstructed theta
#' while the inverse acts on the base-space values.
#' @keywords internal
nsf_apply <- function(net, made, z, x, inverse = FALSE, values = z) {
  p <- net$dim_theta
  pars <- made(z, x)
  n <- z$shape[1]
  flat <- function(t) t$reshape(c(n * p, -1))
  out <- rq_spline(values$reshape(c(n * p)), flat(pars$w), flat(pars$h),
                   flat(pars$d), inverse = inverse,
                   tail_bound = net$tail_bound)
  list(z = out$outputs$reshape(c(n, p)),
       logdet = out$logdet$reshape(c(n, p))$sum(dim = 2))
}

#' Per-row NSF log density (standardized space), as a torch tensor
#' @keywords internal
nsf_log_prob_tensor <- function(net, theta, x) {
  p <- net$dim_theta
  x <- embed_x(net, x)
  rev_idx <- torch::torch_tensor(rev(seq_len(p)), dtype = torch::torch_long())
  z <- theta
  logdet <- torch::torch_zeros(theta$shape[1])
  for (k in seq_len(net$n_transforms)) {
    if (k > 1L && p > 1L) z <- z[, rev_idx, drop = FALSE]
    st <- nsf_apply(net, net$mades[[k]], z, x, inverse = FALSE)
    z <- st$z
    logdet <- logdet + st$logdet
  }
  base_lp <- (-0.5 * (z$pow(2) + log(2 * pi)))$sum(dim = 2)
  base_lp + logdet
}

#' base variable u -> theta, inverting each spline dimension-by-dimension
#' @keywords internal
nsf_inverse <- function(net, u, x) {
  p <- net$dim_theta
  x <- embed_x(net, x)
  rev_idx <- torch::torch_tensor(rev(seq_len(p)), dtype = torch::torch_long())
  z <- u
  for (k in rev(seq_len(net$n_transforms))) {
    made <- net$mades[[k]]
    out <- torch::torch_zeros_like(z)
    for (d in seq_len(p)) {
      st <- nsf_apply(net, made, out, x, inverse = TRUE, values = z)
      out[, d] <- st$z[, d]
    }
    z <- out
    if (k > 1L && p > 1L) z <- z[, rev_idx, drop = FALSE]
  }
  z
}

#' Train an NSF on standardized (theta, x)
#' @keywords internal
fit_nsf <- function(theta, x, n_transforms = 5L, hidden = c(50L, 50L),
                    n_bins = 8L, tail_bound = 3,
                    max_epochs = 500L, batch_size = 100L, lr = 5e-4,
                    validation_fraction = 0.1, patience = 20L,
                    n_restarts = 1L, clip_grad_norm = 5, embedding = NULL,
                    seed = NULL, verbose = FALSE) {
  theta <- as_theta_matrix(theta)
  x <- as_theta_matrix(x)
  dim_theta <- ncol(theta)
  dim_x <- ncol(x)

  trained <- train_conditional_de(
    build_net = function() nsf_module(dim_x, dim_theta, n_transforms, hidden,
                                      n_bins, tail_bound, embedding)(),
    log_prob_fn = nsf_log_prob_tensor,
    theta = theta, x = x,
    max_epochs = max_epochs, batch_size = batch_size, lr = lr,
    validation_fraction = validation_fraction, patience = patience,
    n_restarts = n_restarts, clip_grad_norm = clip_grad_norm,
    seed = seed, verbose = verbose
  )

  structure(
    list(net = trained$net, dim_theta = dim_theta, dim_x = dim_x,
         n_transforms = n_transforms, hidden = hidden, n_bins = n_bins,
         tail_bound = tail_bound, embedding = embedding,
         best_val_loss = trained$best_val_loss, history = trained$history),
    class = c("nsbi_de_nsf", "nsbi_de")
  )
}

#' @export
de_log_prob.nsbi_de_nsf <- function(de, theta, x) {
  theta <- as_theta_matrix(theta, de$dim_theta)
  x <- as_theta_matrix(x, de$dim_x)
  if (nrow(x) == 1L && nrow(theta) > 1L) {
    x <- matrix(x, nrow = nrow(theta), ncol = ncol(x), byrow = TRUE)
  }
  tt <- torch::torch_tensor(theta, dtype = torch::torch_float())
  xt <- torch::torch_tensor(x, dtype = torch::torch_float())
  torch::with_no_grad({
    as.numeric(nsf_log_prob_tensor(de$net, tt, xt)$to(dtype = torch::torch_float64()))
  })
}

#' @export
de_sample.nsbi_de_nsf <- function(de, x, n) {
  x <- as_theta_matrix(x, de$dim_x)[1, , drop = FALSE]
  xrep <- matrix(x, nrow = n, ncol = de$dim_x, byrow = TRUE)
  xt <- torch::torch_tensor(xrep, dtype = torch::torch_float())
  u <- torch::torch_randn(c(n, de$dim_theta))
  torch::with_no_grad({
    torch::as_array(nsf_inverse(de$net, u, xt)$to(dtype = torch::torch_float64()))
  })
}
