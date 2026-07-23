#' Masked Autoregressive Flow (MAF) conditional density estimator
#'
#' A normalizing flow maps parameters \eqn{\theta} to a standard-normal base
#' variable through a stack of invertible transforms, giving exact densities by
#' the change of variables. The MAF (Papamakarios et al., 2017) uses masked
#' autoregressive networks (MADE, Germain et al., 2015): each transform is
#' \deqn{u_d = (\theta_d - \mu_d(\theta_{<d}, x)) \exp(-\alpha_d(\theta_{<d}, x)),}
#' where the masks guarantee that \eqn{\mu_d, \alpha_d} depend only on earlier
#' dimensions of \eqn{\theta} (and freely on the conditioning data `x`). Density
#' evaluation is a single forward pass; sampling inverts the transform one
#' dimension at a time. Between transforms the parameter order is reversed so
#' every dimension gets conditioned on every other across the stack.
#'
#' This is `sbi`'s default flow family and handles non-Gaussian posteriors that
#' the MDN struggles with. Select it with `npe(..., density_estimator = "maf")`.
#'
#' @keywords internal
#' @name maf
NULL

#' Linear layer with a fixed binary mask on the weights.
#' (Defined inside a function so the package loads without torch installed.)
#' @keywords internal
masked_linear <- function(in_features, out_features, mask) {
  torch::nn_module(
    classname = "nsbi_masked_linear",
    initialize = function() {
      self$linear <- torch::nn_linear(in_features, out_features)
      # mask is (out_features, in_features), matching nn_linear's weight
      self$mask <- torch::nn_buffer(
        torch::torch_tensor(mask, dtype = torch::torch_float()))
    },
    forward = function(input) {
      torch::nnf_linear(input, self$linear$weight * self$mask, self$linear$bias)
    }
  )()
}

#' MADE masks for one autoregressive transform.
#'
#' Degrees: theta inputs get 1..p, conditioning inputs get 0 (visible to all),
#' hidden units cycle through 1..(p-1) (or 0 when p = 1). A connection
#' into a hidden unit requires hidden_degree >= input_degree; a connection into
#' output dimension d requires d > hidden_degree. This makes output d a function
#' of theta_{<d} and x only.
#' @keywords internal
made_masks <- function(dim_theta, dim_x, hidden) {
  p <- dim_theta
  in_degrees <- c(seq_len(p), rep(0L, dim_x))
  hidden_choices <- if (p == 1L) 0L else seq_len(p - 1L)
  degrees <- list(in_degrees)
  for (h in hidden) degrees[[length(degrees) + 1L]] <- rep_len(hidden_choices, h)
  masks <- list()
  for (l in seq_along(hidden)) {
    masks[[l]] <- outer(degrees[[l + 1L]], degrees[[l]], `>=`) * 1
  }
  out_mask <- outer(seq_len(p), degrees[[length(degrees)]], `>`) * 1
  list(hidden = masks, out = out_mask)
}

#' One MADE block: (theta, x) -> per-dimension shift mu and log-scale alpha
#' @keywords internal
made_module <- function(dim_x, dim_theta, hidden) {
  torch::nn_module(
    classname = "nsbi_made",
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
      self$out_mu <- masked_linear(sizes[length(sizes)], dim_theta, m$out)
      self$out_alpha <- masked_linear(sizes[length(sizes)], dim_theta, m$out)
      # Start at the identity transform (mu = 0, alpha = 0) for stable training.
      torch::nn_init_zeros_(self$out_mu$linear$weight)
      torch::nn_init_zeros_(self$out_mu$linear$bias)
      torch::nn_init_zeros_(self$out_alpha$linear$weight)
      torch::nn_init_zeros_(self$out_alpha$linear$bias)
    },
    forward = function(theta, x) {
      h <- self$trunk(torch::torch_cat(list(theta, x), dim = 2))
      alpha <- torch::torch_clamp(self$out_alpha(h), min = -8, max = 8)
      list(mu = self$out_mu(h), alpha = alpha)
    }
  )
}

#' The full MAF: a stack of MADE transforms with order reversal in between
#' @keywords internal
maf_module <- function(dim_x, dim_theta, n_transforms, hidden,
                       embedding = NULL) {
  torch::nn_module(
    classname = "nsbi_maf_net",
    initialize = function() {
      self$dim_theta <- dim_theta
      self$n_transforms <- n_transforms
      self$has_embedding <- !is.null(embedding)
      if (self$has_embedding) {
        self$embedding <- build_embedding_module(embedding, dim_x)
      }
      cond_dim <- embedding_output_dim(embedding, dim_x)
      self$mades <- torch::nn_module_list(
        lapply(seq_len(n_transforms),
               function(k) made_module(cond_dim, dim_theta, hidden)())
      )
    },
    forward = function(theta, x) {
      stop("call maf_log_prob_tensor() / maf_sample_tensor() instead")
    }
  )
}

#' theta -> base variable u, accumulating the log |det Jacobian|.
#' Returns list(u = (b, p) tensor, logdet = (b,) tensor).
#' @keywords internal
maf_forward <- function(net, theta, x) {
  p <- net$dim_theta
  x <- embed_x(net, x)
  rev_idx <- torch::torch_tensor(rev(seq_len(p)), dtype = torch::torch_long())
  z <- theta
  logdet <- torch::torch_zeros(theta$shape[1])
  for (k in seq_len(net$n_transforms)) {
    if (k > 1L && p > 1L) z <- z[, rev_idx, drop = FALSE]
    ma <- net$mades[[k]](z, x)
    z <- (z - ma$mu) * torch::torch_exp(-ma$alpha)
    logdet <- logdet - ma$alpha$sum(dim = 2)
  }
  list(u = z, logdet = logdet)
}

#' base variable u -> theta (inverts maf_forward), dimension by dimension
#' @keywords internal
maf_inverse <- function(net, u, x) {
  p <- net$dim_theta
  x <- embed_x(net, x)
  rev_idx <- torch::torch_tensor(rev(seq_len(p)), dtype = torch::torch_long())
  z <- u
  for (k in rev(seq_len(net$n_transforms))) {
    made <- net$mades[[k]]
    out <- torch::torch_zeros_like(z)
    for (d in seq_len(p)) {
      ma <- made(out, x)
      out[, d] <- z[, d] * torch::torch_exp(ma$alpha[, d]) + ma$mu[, d]
    }
    z <- out
    if (k > 1L && p > 1L) z <- z[, rev_idx, drop = FALSE]  # rev inverts itself
  }
  z
}

#' Per-row MAF log density (standardized space), as a torch tensor
#' @keywords internal
maf_log_prob_tensor <- function(net, theta, x) {
  fw <- maf_forward(net, theta, x)
  base_lp <- (-0.5 * (fw$u$pow(2) + log(2 * pi)))$sum(dim = 2)
  base_lp + fw$logdet
}

#' Train a MAF on standardized (theta, x)
#' @keywords internal
fit_maf <- function(theta, x, n_transforms = 5L, hidden = c(50L, 50L),
                    max_epochs = 500L, batch_size = 100L, lr = 5e-4,
                    validation_fraction = 0.1, patience = 20L,
                    n_restarts = 1L, clip_grad_norm = 5, embedding = NULL,
                    seed = NULL, verbose = FALSE) {
  theta <- as_theta_matrix(theta)
  x <- as_theta_matrix(x)
  dim_theta <- ncol(theta)
  dim_x <- ncol(x)

  trained <- train_conditional_de(
    build_net = function()
      maf_module(dim_x, dim_theta, n_transforms, hidden, embedding)(),
    log_prob_fn = maf_log_prob_tensor,
    theta = theta, x = x,
    max_epochs = max_epochs, batch_size = batch_size, lr = lr,
    validation_fraction = validation_fraction, patience = patience,
    n_restarts = n_restarts, clip_grad_norm = clip_grad_norm,
    seed = seed, verbose = verbose
  )

  structure(
    list(net = trained$net, dim_theta = dim_theta, dim_x = dim_x,
         n_transforms = n_transforms, hidden = hidden, embedding = embedding,
         best_val_loss = trained$best_val_loss, history = trained$history),
    class = c("nsbi_de_maf", "nsbi_de")
  )
}

#' @export
de_log_prob.nsbi_de_maf <- function(de, theta, x) {
  theta <- as_theta_matrix(theta, de$dim_theta)
  x <- as_theta_matrix(x, de$dim_x)
  if (nrow(x) == 1L && nrow(theta) > 1L) {
    x <- matrix(x, nrow = nrow(theta), ncol = ncol(x), byrow = TRUE)
  }
  tt <- torch::torch_tensor(theta, dtype = torch::torch_float())
  xt <- torch::torch_tensor(x, dtype = torch::torch_float())
  torch::with_no_grad({
    as.numeric(maf_log_prob_tensor(de$net, tt, xt)$to(dtype = torch::torch_float64()))
  })
}

#' @export
de_sample.nsbi_de_maf <- function(de, x, n) {
  x <- as_theta_matrix(x, de$dim_x)[1, , drop = FALSE]
  xrep <- matrix(x, nrow = n, ncol = de$dim_x, byrow = TRUE)
  xt <- torch::torch_tensor(xrep, dtype = torch::torch_float())
  u <- torch::torch_randn(c(n, de$dim_theta))
  torch::with_no_grad({
    torch::as_array(maf_inverse(de$net, u, xt)$to(dtype = torch::torch_float64()))
  })
}
