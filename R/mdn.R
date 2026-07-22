#' Mixture Density Network (MDN) conditional density estimator
#'
#' The MDN is the default neural density estimator in `neuralsbi`. A multilayer
#' perceptron maps the data `x` to the parameters of a Gaussian mixture over the
#' parameters \eqn{\theta}: mixture logits, component means, and (full)
#' lower-triangular Cholesky factors of each component covariance. Training
#' minimizes the negative log-likelihood of \eqn{\theta} under the mixture,
#' which -- when simulations are drawn from the prior -- yields a direct
#' amortized approximation of the posterior \eqn{p(\theta \mid x)}.
#'
#' A native R/`torch` implementation of the multivariate-Gaussian mixture
#' density network (Bishop, 1994).
#'
#' @keywords internal
#' @name mdn
NULL

#' @keywords internal
tril_indices <- function(p) {
  idx <- which(lower.tri(matrix(0, p, p), diag = TRUE), arr.ind = TRUE)
  idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]  # row-major
  list(row = idx[, 1], col = idx[, 2], is_diag = idx[, 1] == idx[, 2])
}

#' Build the MDN torch module
#' @keywords internal
mdn_module <- function(dim_x, dim_theta, n_components, hidden,
                       embedding = NULL) {
  torch::nn_module(
    classname = "nsbi_mdn_net",
    initialize = function() {
      self$dim_theta <- dim_theta
      self$K <- n_components
      self$tri <- tril_indices(dim_theta)
      self$tril_size <- length(self$tri$row)
      self$has_embedding <- !is.null(embedding)
      if (self$has_embedding) {
        self$embedding <- build_embedding_module(embedding, dim_x)
      }
      layers <- list()
      prev <- embedding_output_dim(embedding, dim_x)
      for (h in hidden) {
        layers[[length(layers) + 1L]] <- torch::nn_linear(prev, h)
        layers[[length(layers) + 1L]] <- torch::nn_relu()
        prev <- h
      }
      self$trunk <- do.call(torch::nn_sequential, layers)
      self$head_logits <- torch::nn_linear(prev, self$K)
      self$head_means <- torch::nn_linear(prev, self$K * dim_theta)
      self$head_tril <- torch::nn_linear(prev, self$K * self$tril_size)
    },
    forward = function(x) {
      if (self$has_embedding) x <- self$embedding(x)
      h <- self$trunk(x)
      batch <- x$shape[1]
      logits <- self$head_logits(h)                                # (b, K)
      means <- self$head_means(h)$view(c(batch, self$K, self$dim_theta))
      tril_flat <- self$head_tril(h)$view(c(batch, self$K, self$tril_size))
      list(logits = logits, means = means, tril_flat = tril_flat)
    }
  )
}

#' Assemble batched lower-triangular Cholesky factors from the flat head output.
#' Diagonal entries are passed through softplus (+ eps) to stay positive.
#' Returns a tensor of shape (batch, K, p, p).
#' @keywords internal
mdn_build_tril <- function(net, tril_flat) {
  b <- tril_flat$shape[1]
  K <- net$K
  p <- net$dim_theta
  L <- torch::torch_zeros(c(b, K, p, p), dtype = tril_flat$dtype)
  for (m in seq_len(net$tril_size)) {
    val <- tril_flat[, , m]
    if (net$tri$is_diag[m]) {
      val <- torch::nnf_softplus(val) + 1e-6
    }
    L[, , net$tri$row[m], net$tri$col[m]] <- val
  }
  L
}

#' Per-row mixture log density (in standardized theta space), as a torch tensor.
#' `theta`: (b, p) tensor, `x`: (b, q) tensor.
#' @keywords internal
mdn_log_prob_tensor <- function(net, theta, x) {
  params <- net(x)
  L <- mdn_build_tril(net, params$tril_flat)                    # (b,K,p,p)
  p <- net$dim_theta
  diff <- theta$unsqueeze(2) - params$means                     # (b,K,p)
  diff_col <- diff$unsqueeze(4)                                  # (b,K,p,1)
  z <- torch::linalg_solve_triangular(L, diff_col, upper = FALSE) # (b,K,p,1)
  quad <- z$pow(2)$sum(dim = c(3, 4))                            # (b,K)
  diag_L <- torch::torch_diagonal(L, dim1 = 3, dim2 = 4)        # (b,K,p)
  logdet <- 2 * torch::torch_log(diag_L)$sum(dim = 3)           # (b,K)
  const <- p * log(2 * pi)
  comp_lp <- -0.5 * (const + logdet + quad)                     # (b,K)
  log_w <- torch::nnf_log_softmax(params$logits, dim = 2)       # (b,K)
  torch::torch_logsumexp(log_w + comp_lp, dim = 2)              # (b,)
}

#' Train an MDN on standardized (theta, x)
#' @keywords internal
fit_mdn <- function(theta, x, n_components = 5L, hidden = c(50L, 50L),
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
      mdn_module(dim_x, dim_theta, n_components, hidden, embedding)(),
    log_prob_fn = mdn_log_prob_tensor,
    theta = theta, x = x,
    max_epochs = max_epochs, batch_size = batch_size, lr = lr,
    validation_fraction = validation_fraction, patience = patience,
    n_restarts = n_restarts, clip_grad_norm = clip_grad_norm,
    seed = seed, verbose = verbose
  )

  structure(
    list(net = trained$net, dim_theta = dim_theta, dim_x = dim_x,
         n_components = n_components, hidden = hidden, embedding = embedding,
         best_val_loss = trained$best_val_loss, history = trained$history),
    class = c("nsbi_de_mdn", "nsbi_de")
  )
}

#' @export
de_log_prob.nsbi_de_mdn <- function(de, theta, x) {
  theta <- as_theta_matrix(theta, de$dim_theta)
  x <- as_theta_matrix(x, de$dim_x)
  if (nrow(x) == 1L && nrow(theta) > 1L) {
    x <- matrix(x, nrow = nrow(theta), ncol = ncol(x), byrow = TRUE)
  }
  tt <- torch::torch_tensor(theta, dtype = torch::torch_float())
  xt <- torch::torch_tensor(x, dtype = torch::torch_float())
  torch::with_no_grad({
    as.numeric(mdn_log_prob_tensor(de$net, tt, xt)$to(dtype = torch::torch_float64()))
  })
}

#' @export
de_sample.nsbi_de_mdn <- function(de, x, n) {
  x <- as_theta_matrix(x, de$dim_x)[1, , drop = FALSE]
  xt <- torch::torch_tensor(x, dtype = torch::torch_float())
  params <- torch::with_no_grad(de$net(xt))
  L <- torch::with_no_grad(mdn_build_tril(de$net, params$tril_flat))
  logits <- as.numeric(torch::as_array(params$logits[1, ]))
  means <- torch::as_array(params$means[1, , ])           # K x p
  Larr <- torch::as_array(L[1, , , ])                      # K x p x p
  if (de$n_components == 1L) {
    means <- matrix(means, nrow = 1L)
    Larr <- array(Larr, dim = c(1L, de$dim_theta, de$dim_theta))
  }
  w <- exp(logits - max(logits)); w <- w / sum(w)
  comp <- sample.int(de$n_components, n, replace = TRUE, prob = w)
  out <- matrix(0, nrow = n, ncol = de$dim_theta)
  z <- matrix(stats::rnorm(n * de$dim_theta), nrow = n)
  for (k in unique(comp)) {
    rows <- which(comp == k)
    Lk <- matrix(Larr[k, , ], de$dim_theta, de$dim_theta)
    out[rows, ] <- sweep(z[rows, , drop = FALSE] %*% t(Lk), 2, means[k, ], `+`)
  }
  out
}
