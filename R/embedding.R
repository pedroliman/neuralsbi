#' Embedding (summary) networks for structured observations
#'
#' Raw observations are often high-dimensional or structured (a time series, a
#' set of summary statistics, an image) where feeding `x` straight into the
#' density estimator wastes capacity. An embedding network learns a low-
#' dimensional summary \eqn{h = f_\psi(x)} jointly with the density estimator,
#' so the conditioning path becomes \eqn{q_\phi(\theta \mid f_\psi(x))}. This
#' mirrors `sbi`'s `embedding_net` argument.
#'
#' `embedding_mlp()` builds a multilayer-perceptron summary network: a stack of
#' fully connected ReLU layers mapping the (standardized) data to a vector of
#' `output_dim` features. Pass the result to [npe()] via `embedding_net`; it is
#' trained end to end with the estimator and its parameters live inside the
#' fitted network, so sampling and `log_prob` route through it automatically.
#'
#' The embedding consumes the standardized data (the same z-scoring [npe()]
#' applies to `x` without an embedding), which keeps the summary network's
#' inputs on a common scale. Standardization of the *features* is intentionally
#' left to the network itself; the estimators operate on the raw embedding
#' output.
#'
#' @param output_dim Number of summary features the network emits. This is the
#'   effective data dimension the density estimator conditions on.
#' @param hidden Integer vector of hidden-layer widths (ReLU between layers).
#'   An empty vector gives a single linear map to `output_dim`.
#'
#' @return An `nsbi_embedding` specification. It carries no torch objects (the
#'   network is built lazily at fit time), so it is safe to construct without
#'   `torch` installed.
#'
#' @examples
#' emb <- embedding_mlp(output_dim = 8, hidden = c(64, 64))
#' # fit <- npe(prior, simulator, density_estimator = "maf", embedding_net = emb)
#' @seealso [npe()]
#' @export
embedding_mlp <- function(output_dim = 16L, hidden = c(64L, 64L)) {
  if (length(output_dim) != 1L || !is.finite(output_dim) || output_dim < 1L) {
    stop("`output_dim` must be a single positive integer.", call. = FALSE)
  }
  hidden <- as.integer(hidden)
  if (length(hidden) && any(!is.finite(hidden) | hidden < 1L)) {
    stop("`hidden` must be positive integers.", call. = FALSE)
  }
  structure(
    list(type = "mlp", output_dim = as.integer(output_dim), hidden = hidden),
    class = "nsbi_embedding"
  )
}

#' Effective conditioning dimension after an (optional) embedding.
#'
#' The identity embedding (`spec = NULL`) leaves the data dimension unchanged;
#' otherwise the estimator conditions on `output_dim` features.
#' @keywords internal
embedding_output_dim <- function(spec, dim_x) {
  if (is.null(spec)) dim_x else spec$output_dim
}

#' Build the embedding torch submodule, or `NULL` for the identity embedding.
#'
#' Constructed lazily so no torch object exists at package-load time. Returns an
#' instantiated `nn_module` (call site stores it as a submodule so its
#' parameters train jointly and travel with the estimator's `state_dict`).
#' @keywords internal
build_embedding_module <- function(spec, dim_x) {
  if (is.null(spec)) return(NULL)
  hidden <- spec$hidden
  output_dim <- spec$output_dim
  torch::nn_module(
    classname = "nsbi_embedding_mlp",
    initialize = function() {
      layers <- list()
      prev <- dim_x
      for (h in hidden) {
        layers[[length(layers) + 1L]] <- torch::nn_linear(prev, h)
        layers[[length(layers) + 1L]] <- torch::nn_relu()
        prev <- h
      }
      layers[[length(layers) + 1L]] <- torch::nn_linear(prev, output_dim)
      self$net <- do.call(torch::nn_sequential, layers)
    },
    forward = function(x) self$net(x)
  )()
}

#' Apply an estimator's embedding to conditioning data, if it has one.
#'
#' Estimators call this once per forward/inverse pass so the embedding runs a
#' single time; the raw (standardized) `x` still enters at the `de_*` boundary,
#' keeping `de$dim_x` the raw data dimension.
#' @keywords internal
embed_x <- function(net, x) {
  if (isTRUE(net$has_embedding)) net$embedding(x) else x
}
