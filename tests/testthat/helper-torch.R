# Shared across test files: skip neural tests when libtorch is unavailable.
skip_if_no_torch <- function() {
  testthat::skip_if_not(
    requireNamespace("torch", quietly = TRUE) &&
      isTRUE(torch::torch_is_installed()),
    "torch not available"
  )
}
