#' @keywords internal
"_PACKAGE"

# `self` is injected by torch::nn_module() inside initialize()/forward();
# declare it to silence a spurious "no visible binding" NOTE.
utils::globalVariables("self")

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

# Cache for the libtorch load probe. Loading the native library is expensive
# and side-effecting, so probe at most once per session and remember both the
# outcome and, on failure, the underlying error for a targeted message.
.torch_load_state <- new.env(parent = emptyenv())

# Try to actually load the libtorch native library.
#
# `torch::torch_is_installed()` only checks that the library *files* are on
# disk; it never dlopen()s them. On macOS the bundled libtorch is sometimes
# built against a newer macOS than the one running (an "image not found" /
# "Symbol not found" dlopen failure), so the files are present but any tensor
# op crashes with a cryptic "Lantern is not loaded". Force the load once by
# touching a tensor and cache the result.
#' @keywords internal
torch_loadable <- function() {
  if (!requireNamespace("torch", quietly = TRUE)) return(FALSE)
  if (!isTRUE(torch::torch_is_installed())) return(FALSE)
  if (!is.null(.torch_load_state$ok)) return(.torch_load_state$ok)
  ok <- tryCatch(
    suppressWarnings(suppressMessages({
      torch::torch_tensor(1)
      TRUE
    })),
    error = function(e) {
      .torch_load_state$error <- conditionMessage(e)
      FALSE
    }
  )
  .torch_load_state$ok <- ok
  ok
}

# Build a clear, platform-aware message for an installed-but-unloadable torch.
#' @keywords internal
torch_load_error <- function() {
  detail <- .torch_load_state$error %||%
    "the libtorch native library could not be loaded."
  lines <- c(
    "'torch' is installed but its native library (libtorch) failed to load.",
    paste0("  Underlying error: ", detail),
    ""
  )
  if (identical(Sys.info()[["sysname"]], "Darwin")) {
    lines <- c(
      lines,
      "On macOS this usually means the bundled libtorch was built for a newer",
      "macOS than the one you are running (note the 'newer than running OS' or",
      "'Symbol not found' text above). Reinstalling the same torch will not",
      "help. Instead, either:",
      "  1. Update macOS to the version libtorch was built for, or",
      "  2. Install an earlier 'torch' whose libtorch targets your macOS, e.g.:",
      "       remotes::install_version('torch', '0.13.0')",
      "       torch::install_torch(reinstall = TRUE)",
      "     then restart R. Try progressively older versions until it loads.",
      "See the macOS troubleshooting notes:",
      "  https://github.com/pedroliman/neuralsbi#troubleshooting-torch-on-macos"
    )
  } else {
    lines <- c(
      lines,
      "Try reinstalling libtorch and restarting R:",
      "  torch::install_torch(reinstall = TRUE)"
    )
  }
  paste(lines, collapse = "\n")
}

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
  if (!torch_loadable()) {
    stop(torch_load_error(), call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
torch_available <- function() {
  requireNamespace("torch", quietly = TRUE) &&
    isTRUE(torch::torch_is_installed()) &&
    torch_loadable()
}

#' @keywords internal
verbose_cat <- function(verbose, ...) {
  if (isTRUE(verbose)) cat(...)
}
