# Shared across test files: skip neural tests when libtorch is unavailable.
skip_if_no_torch <- function() {
  testthat::skip_if_not(neuralsbi:::torch_available(), "torch not available")
}

# c2st()'s default classifier is sbibm's MLP, which trains on torch. The
# analytic-parity tests still want a C2ST on a machine without libtorch, so
# they assert the torch-free logistic test always and add the MLP where it can
# run. Callers use this rather than skipping the whole test.
has_torch <- function() isTRUE(neuralsbi:::torch_available())
