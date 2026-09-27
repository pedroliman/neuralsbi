# Shared plumbing for the sbibm benchmark scripts: argument parsing, paths,
# logging, and loading neuralsbi itself.

bench_dir <- function() {
  # Works both when sourced from a script run with Rscript and interactively
  # from the package root.
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg))))
  }
  if (dir.exists("dev/benchmarks")) return(normalizePath("dev/benchmarks"))
  normalizePath(".")
}

BENCH_DIR <- bench_dir()

pkg_root <- function() normalizePath(file.path(BENCH_DIR, "..", ".."))

#' Load neuralsbi, preferring the installed package and falling back to the
#' working tree so the scripts also run against uncommitted changes.
load_neuralsbi <- function() {
  if (requireNamespace("neuralsbi", quietly = TRUE)) {
    suppressPackageStartupMessages(library(neuralsbi))
    return(invisible("installed"))
  }
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(pkg_root(), quiet = TRUE)
    return(invisible("load_all"))
  }
  stop("neuralsbi is not installed and pkgload is unavailable. Run\n",
       "  R CMD INSTALL ", pkg_root(), call. = FALSE)
}

#' Minimal `--key value` / `--key=value` command line parser.
#'
#' @param defaults Named list of defaults. The type of each default decides how
#'   the supplied value is coerced; a comma-separated value becomes a vector.
parse_args <- function(defaults, args = commandArgs(trailingOnly = TRUE)) {
  out <- defaults
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (!startsWith(a, "--")) {
      stop("Unexpected argument: ", a, call. = FALSE)
    }
    if (grepl("=", a, fixed = TRUE)) {
      key <- sub("^--", "", sub("=.*$", "", a))
      val <- sub("^[^=]*=", "", a)
    } else {
      key <- sub("^--", "", a)
      nxt <- if (i < length(args)) args[[i + 1]] else NA_character_
      if (is.na(nxt) || startsWith(nxt, "--")) {
        val <- "TRUE"
      } else {
        val <- nxt
        i <- i + 1
      }
    }
    key <- gsub("-", "_", key)
    if (!key %in% names(defaults)) {
      stop("Unknown option --", key, ". Known options: ",
           paste(names(defaults), collapse = ", "), call. = FALSE)
    }
    out[[key]] <- coerce_like(val, defaults[[key]])
    i <- i + 1
  }
  out
}

coerce_like <- function(val, template) {
  parts <- strsplit(val, ",", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  if (is.logical(template)) return(as.logical(parts))
  if (is.numeric(template)) return(as.numeric(parts))
  parts
}

say <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
  utils::flush.console()
}

#' Format a simulation budget the way the paper's tables do.
budget_label <- function(n) {
  e <- round(log10(n))
  if (10^e == n) sprintf("10^%d", e) else format(n, scientific = FALSE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
