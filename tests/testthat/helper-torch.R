# Shared across test files: skip neural tests when libtorch is unavailable.
# "Unavailable" covers both a missing install and one that is present on disk
# but fails to load (e.g. a libtorch built for a newer macOS), so a broken
# back end skips cleanly instead of erroring the suite.
skip_if_no_torch <- function() {
  available <- requireNamespace("torch", quietly = TRUE) &&
    isTRUE(torch::torch_is_installed())
  # torch_loadable() is internal to neuralsbi; it is visible when tests run in
  # the package namespace (R CMD check, devtools::test()).
  if (available && exists("torch_loadable", mode = "function")) {
    available <- isTRUE(torch_loadable())
  }
  testthat::skip_if_not(available, "torch not available")
}
