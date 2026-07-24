#' @keywords internal
"_PACKAGE"

# `self` is injected by torch::nn_module() inside initialize()/forward();
# `.data` is rlang/ggplot2's tidy-eval pronoun, used inside aes() in
# plotting.R. Declare both to silence spurious "no visible binding" NOTEs.
# ggplot2 is a Suggests, not an Imports, so this -- not @importFrom -- is the
# dependency-safe way to quiet the check.
utils::globalVariables(c("self", ".data"))

#' Coerce parameters/data to a numeric matrix with a known column count
#' @keywords internal
as_theta_matrix <- function(x, d = NULL) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x))) {
    # a plain vector: interpret as a single row if length matches d,
    # otherwise as a column of 1-D values.
    if (!is.null(d) && length(x) == d) {
      x <- matrix(x, nrow = 1L)
    } else {
      x <- matrix(x, ncol = if (is.null(d)) 1L else d, byrow = TRUE)
    }
  }
  storage.mode(x) <- "double"
  if (!is.null(d) && ncol(x) != d) {
    stop(sprintf("Expected %d columns but got %d.", d, ncol(x)), call. = FALSE)
  }
  x
}

#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Check that torch is available, error otherwise
#' @keywords internal
require_torch <- function() {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop(
      "This density estimator needs the 'torch' package.\n",
      "Install it with install.packages('torch') and then torch::install_torch().\n",
      "Alternatively use density_estimator = 'linear_gaussian' for a torch-free baseline.",
      call. = FALSE
    )
  }
  if (!torch::torch_is_installed()) {
    stop(
      "'torch' is installed but libtorch is not. Run torch::install_torch().",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @keywords internal
torch_available <- function() {
  requireNamespace("torch", quietly = TRUE) && isTRUE(torch::torch_is_installed())
}

#' Check that ggplot2 (and, for [pairplot()], GGally) are available
#' @keywords internal
require_ggplot2 <- function(ggally = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "This plotting function needs the 'ggplot2' package.\n",
      "Install it with install.packages('ggplot2').",
      call. = FALSE
    )
  }
  if (ggally && !requireNamespace("GGally", quietly = TRUE)) {
    stop(
      "pairplot() needs the 'GGally' package (for ggpairs()).\n",
      "Install it with install.packages('GGally').",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @keywords internal
verbose_cat <- function(verbose, ...) {
  if (isTRUE(verbose)) cat(...)
}
