# Shared across test files: skip neural tests when libtorch is unavailable.
skip_if_no_torch <- function() {
  testthat::skip_if_not(neuralsbi:::torch_available(), "torch not available")
}
