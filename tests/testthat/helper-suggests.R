# Shared across test files: skip plotting tests when the graphics Suggests are
# absent. ggplot2, GGally and ggdensity are Suggests, so `R CMD check` with
# `_R_CHECK_FORCE_SUGGESTS_=false` -- the way CRAN checks a package on a
# machine that has not installed them -- must not fail here. Same contract as
# skip_if_no_torch() in helper-torch.R: the suite runs everywhere.
skip_if_no_ggplot2 <- function() {
  testthat::skip_if_not_installed("ggplot2")
}

skip_if_no_ggally <- function() {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("GGally")
  testthat::skip_if_not_installed("ggdensity")
}

# posterior is a Suggests, used only to cross-check mcmc_diagnostics() against
# a reference implementation. Same contract as skip_if_no_torch().
skip_if_no_posterior <- function() {
  testthat::skip_if_not_installed("posterior")
}
