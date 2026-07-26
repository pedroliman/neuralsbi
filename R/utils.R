#' @keywords internal
"_PACKAGE"

# `self` is injected by torch::nn_module() inside initialize()/forward();
# `.data` is rlang/ggplot2's tidy-eval pronoun, used inside aes() in
# plotting.R. Declare both to silence spurious "no visible binding" NOTEs.
# ggplot2 is a Suggests, not an Imports, so this -- not @importFrom -- is the
# dependency-safe way to quiet the check.
utils::globalVariables(c("self", ".data"))

#' Coerce parameters/data to a numeric matrix with a known column count
#'
#' Preserves column names where they carry meaning: a data frame's or matrix's
#' existing `colnames`, or a plain vector's `names()` when the vector is
#' interpreted as a single row. A vector re-interpreted as a stacked column
#' (the `byrow` branch) loses its names -- they described entries, not a
#' shared parameter/outcome identity.
#' @keywords internal
as_theta_matrix <- function(x, d = NULL) {
  nm <- if (is.data.frame(x) || !is.null(dim(x))) colnames(x) else names(x)
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x))) {
    # a plain vector: interpret as a single row if length matches d,
    # otherwise as a column of 1-D values.
    if (!is.null(d) && length(x) == d) {
      x <- matrix(x, nrow = 1L)
    } else {
      x <- matrix(x, ncol = if (is.null(d)) 1L else d, byrow = TRUE)
      nm <- NULL
    }
  }
  storage.mode(x) <- "double"
  if (!is.null(d) && ncol(x) != d) {
    stop(sprintf("Expected %d columns but got %d.", d, ncol(x)), call. = FALSE)
  }
  if (is.null(colnames(x)) && !is.null(nm) && length(nm) == ncol(x)) {
    colnames(x) <- nm
  }
  x
}

#' Parse a label as a plotmath expression when possible
#'
#' Parameter/outcome names that happen to be valid R syntax (`"beta[1]"`,
#' `"rho"`, `"sigma^2"`) render as their mathematical symbol -- Greek letters,
#' sub/superscripts -- when passed through R's plotmath. Names that are not
#' parseable, or `NULL`, pass through unchanged so callers can fall back to a
#' plain-text label.
#' @keywords internal
math_expr <- function(label) {
  if (is.null(label)) return(NULL)
  expr <- tryCatch(str2lang(label), error = function(e) NULL)
  expr %||% label
}

#' Vectorized, parse-safe plotmath label text
#'
#' Re-quotes any entry that is not valid R syntax as a string literal, which
#' always parses (and renders as plain text under plotmath). Applying this
#' before a label is used as a facet/column name guarantees
#' `ggplot2::label_parsed()` never hits a parse error, whether or not the
#' original label happens to look like math (`"beta[1]"`) or not (`"growth
#' rate"`).
#' @keywords internal
math_safe_text <- function(labels) {
  labels <- as.character(labels)
  vapply(labels, function(s) {
    ok <- tryCatch({ parse(text = s); TRUE }, error = function(e) FALSE)
    if (ok) s else paste0('"', gsub('"', '\\"', s, fixed = TRUE), '"')
  }, character(1), USE.NAMES = FALSE)
}

#' Vectorized, parse-safe plotmath labels
#'
#' Like [math_expr()] but for a whole vector at once, returning an
#' `expression()` (so a discrete ggplot2 scale's `labels =` can render each
#' entry as plotmath, mirroring `ggplot2::label_parsed()` for facet strips).
#' @keywords internal
math_labels <- function(labels) {
  parse(text = math_safe_text(labels))
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

#' Check that ggplot2 (and, for [pairplot()], GGally and ggdensity) are available
#' @keywords internal
require_ggplot2 <- function(ggally = FALSE, ggdensity = FALSE) {
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
  if (ggdensity && !requireNamespace("ggdensity", quietly = TRUE)) {
    stop(
      "pairplot() needs the 'ggdensity' package (for geom_hdr()).\n",
      "Install it with install.packages('ggdensity').",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @keywords internal
verbose_cat <- function(verbose, ...) {
  if (isTRUE(verbose)) cat(...)
}
