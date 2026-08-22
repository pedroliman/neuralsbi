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
#'
#' @param what What needs torch, as the subject of the sentence. The default
#'   suits the density estimators, which is where most of the calls are.
#' @param alternative The torch-free thing to do instead, as a full sentence.
#'   Every caller has one, and naming it is the difference between a dead end
#'   and a next step.
#' @keywords internal
require_torch <- function(what = "This density estimator",
                          alternative = paste("Alternatively use",
                                              "density_estimator =",
                                              "'linear_gaussian' for a",
                                              "torch-free baseline.")) {
  if (torch_available()) return(invisible(TRUE))
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop(
      what, " needs the 'torch' package.\n",
      "Install it with install.packages('torch') and then torch::install_torch().\n",
      alternative,
      call. = FALSE
    )
  }
  stop(
    "'torch' is installed but libtorch is not. Run torch::install_torch().",
    call. = FALSE
  )
}

#' @keywords internal
torch_available <- function() {
  requireNamespace("torch", quietly = TRUE) && isTRUE(torch::torch_is_installed())
}

#' Resolve a validated `device` keyword to a concrete, available torch device
#'
#' [check_device_arg()] has already rejected anything but `"cpu"`, `"cuda"`,
#' `"mps"`, `"gpu"` or `"auto"`; this turns that keyword into `"cpu"`,
#' `"cuda"` or `"mps"`, matching what Python `sbi`'s `process_device()` does.
#' Requires `torch` (call [require_torch()] first; every caller here is about
#' to build a network, so it already has).
#'
#' `"cuda"`/`"mps"` name a specific device, so asking for one that is not
#' there errors rather than downgrading silently -- a silent fallback would
#' hide a real problem (a missing CUDA build, a non-Apple-silicon Mac).
#' `"gpu"`/`"auto"` never named a specific device, so it resolves
#' CUDA -> MPS -> CPU and *can* fall back silently, mirroring `sbi`'s `"gpu"`.
#'
#' @param device One of `"cpu"`, `"cuda"`, `"mps"`, `"gpu"`, `"auto"`.
#' @return `"cpu"`, `"cuda"` or `"mps"`.
#' @keywords internal
resolve_device <- function(device) {
  if (device == "cpu") return("cpu")
  if (device %in% c("gpu", "auto")) {
    if (isTRUE(torch::cuda_is_available())) return("cuda")
    if (isTRUE(torch::backends_mps_is_available())) return("mps")
    return("cpu")
  }
  if (device == "cuda" && !isTRUE(torch::cuda_is_available())) {
    stop("`device = \"cuda\"` was requested, but torch::cuda_is_available() ",
         "is FALSE on this machine.\nUse device = \"cpu\" (the default), or ",
         "\"gpu\"/\"auto\" to fall back to whichever device is available.",
         call. = FALSE)
  }
  if (device == "mps" && !isTRUE(torch::backends_mps_is_available())) {
    stop("`device = \"mps\"` was requested, but ",
         "torch::backends_mps_is_available() is FALSE on this machine.\n",
         "Use device = \"cpu\" (the default), or \"gpu\"/\"auto\" to fall ",
         "back to whichever device is available.",
         call. = FALSE)
  }
  device
}

#' The torch device a fitted net's parameters currently live on
#'
#' `de_log_prob()`/`de_sample()` need to know where to put their input
#' tensors, and the one place that is always true is the net itself --
#' `de$device` (set at fit time) is a record of that, not a guarantee, since
#' nothing stops a user from moving the net with `net$to()` afterwards. A net
#' with no parameters (there is no such estimator today, but nothing rules one
#' out) defaults to CPU.
#' @keywords internal
net_device <- function(net) {
  params <- net$parameters
  if (length(params)) params[[1L]]$device else torch::torch_device("cpu")
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

#' Print the fit-summary block shared by `print.nsbi_npe()`, `print.nsbi_nle()`
#' and `print.nsbi_snpe()`
#'
#' Parameter names, outcome names, the embedding line (when the fit used
#' one), the simulation count and any drops, the best validation loss, and
#' the dead-network warning read the same fields regardless of which
#' factorization was learned. `save_fn_name` is the one line that legitimately
#' differs between an NPE and an NLE fit, so it is the one argument callers
#' must supply.
#'
#' @param x An `nsbi_npe`- or `nsbi_nle`-family fit.
#' @param save_fn_name Name of the save function to point to in the
#'   dead-network warning, e.g. `"save_npe"`; the matching `load_*()` name is
#'   derived from it.
#' @param data_suffix Text appended to the "data (dim)" line before its
#'   newline, e.g. `"  per observation"` for [nle()].
#' @keywords internal
cat_fit_common <- function(x, save_fn_name, data_suffix = "") {
  cat(sprintf("  parameters (dim)  : %d\n", x$dim_theta))
  if (!is.null(x$param_names)) {
    cat("    names           :", paste(x$param_names, collapse = ", "), "\n")
  }
  cat(sprintf("  data (dim)        : %d%s\n", x$dim_x, data_suffix))
  if (!is.null(x$x_names)) {
    cat("    names           :", paste(x$x_names, collapse = ", "), "\n")
  }
  if (!is.null(x$de$embedding)) {
    cat(sprintf("  embedding (mlp)   : %d -> %d features\n",
                x$dim_x, x$de$embedding$output_dim))
  }
  cat(sprintf("  simulations       : %d\n", x$n_simulations))
  cat_dropped(x$n_dropped, x$n_simulations + x$n_dropped, "non-finite output")
  if (!is.null(x$de$best_val_loss) && is.finite(x$de$best_val_loss)) {
    cat(sprintf("  best val loss     : %.4f\n", x$de$best_val_loss))
  }
  if (!torch_net_alive(x$de$net)) {
    cat("  ! network unusable: a torch fit does not survive saveRDS();\n")
    cat(sprintf("    save with %s() and reload with %s().\n",
                save_fn_name, sub("^save_", "load_", save_fn_name)))
  }
}

#' Print an "N dropped" line shared by fit summaries and calibration
#' diagnostics
#'
#' A fit knows its simulation budget, so its line reports a drop rate;
#' [sbc()] and [tarp()] only know how many further trials were lost mid-run,
#' so leaving `total` at `NULL` switches to that shorter wording.
#'
#' @param n_dropped Number dropped. Nothing is printed when this is `NULL`
#'   or `0`.
#' @param total Simulations attempted (dropped plus kept), or `NULL` for the
#'   diagnostics wording, which has no rate to report.
#' @param what What was dropped, e.g. `"non-finite output"`.
#' @keywords internal
cat_dropped <- function(n_dropped, total = NULL, what) {
  if (is.null(n_dropped) || n_dropped == 0L) return(invisible())
  if (is.null(total)) {
    cat(sprintf("  %d %s\n", n_dropped, what))
  } else {
    cat(sprintf("    dropped         : %d of %d, %s (%.1f%%)\n",
                n_dropped, total, what, 100 * n_dropped / total))
  }
  invisible()
}
