# Shared across test files: skip plotting tests when the graphics Suggests are
# absent. ggplot2 and GGally are Suggests, so `R CMD check` with
# `_R_CHECK_FORCE_SUGGESTS_=false` -- the way CRAN checks a package on a
# machine that has not installed them -- must not fail here. Same contract as
# skip_if_no_torch() in helper-torch.R: the suite runs everywhere.
skip_if_no_ggplot2 <- function() {
  testthat::skip_if_not_installed("ggplot2")
}

skip_if_no_ggally <- function() {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("GGally")
}
