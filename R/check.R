#' Argument validation shared across the package
#'
#' One place for the checks that public entry points repeat. Every message
#' names the argument it is about, says what was actually wrong, and where
#' there is a relevant help topic it points at it. That is the voice of
#' [as_sim_draw()], which is the model these follow: a user who mistypes an
#' argument should not have to guess which of `theta` and `x` the complaint is
#' about, or read a coercion error raised three frames down.
#'
#' These are internal helpers, not exported. Call them at public boundaries;
#' internal code that already knows its shapes keeps using [as_theta_matrix()].
#'
#' @name check
#' @keywords internal
NULL

#' Describe a value in an error message
#'
#' Short and literal: the value itself when it is one number, and otherwise
#' what it is instead of a number.
#' @param x The value to describe.
#' @keywords internal
describe_value <- function(x) {
  if (is.null(x)) return("NULL")
  if (length(x) != 1L) {
    return(sprintf("a length-%d %s vector", length(x), class(x)[1L]))
  }
  if (!is.numeric(x) && !is.logical(x)) {
    return(sprintf("a %s value", class(x)[1L]))
  }
  if (is.nan(x)) return("NaN")
  if (is.na(x)) return("NA")
  format(x)
}

#' Count a noun, pluralizing it
#' @param n The count.
#' @param word The singular noun.
#' @keywords internal
n_things <- function(n, word) {
  sprintf("%d %s%s", n, word, if (n == 1L) "" else "s")
}

#' Describe one parameter set in an error message
#'
#' The values that produced a failure, as `mu = 1.83, sigma = 0.42`, or bare
#' numbers when the prior does not name its parameters. Rounded, because four
#' significant digits are enough to recognise a draw and a full-precision
#' double is not. Only the first `max_show` are printed so a 40-parameter model
#' does not fill the console.
#'
#' @param theta_i One parameter set, named or not.
#' @param max_show How many values to print before truncating.
#' @keywords internal
describe_params <- function(theta_i, max_show = 6L) {
  n <- length(theta_i)
  if (n == 0L) return("none")
  keep <- seq_len(min(n, max_show))
  shown <- vapply(theta_i[keep], function(z) {
    describe_value(if (is.numeric(z)) signif(z, 4L) else z)
  }, character(1))
  nm <- names(theta_i)
  if (!is.null(nm) && all(nzchar(nm[keep]))) {
    shown <- paste(nm[keep], shown, sep = " = ")
  }
  sprintf("%s%s", paste(unname(shown), collapse = ", "),
          if (n > max_show) sprintf(", ... (%s in all)",
                                    n_things(n, "parameter")) else "")
}

#' Require numeric data, naming any column that is not
#'
#' The type half of [check_matrix()], split out because the shape rules differ
#' between entry points but this rule does not. On the pre-computed
#' `(theta, x)` path a row is one simulation, so a bare vector is a column of
#' values rather than [check_matrix()]'s single row; a character or factor
#' column is the same mistake in either place. Left to
#' `storage.mode(x) <- "double"` such a column becomes all `NA`, every row is
#' then dropped as non-finite, and the error blames a simulator that was never
#' called.
#'
#' @param value The user's value: a numeric vector, matrix or data frame.
#' @param arg Name of the argument, as it appears in the user's call.
#' @return `value`, with a data frame converted to a matrix.
#' @keywords internal
check_numeric <- function(value, arg) {
  bad <- function(fmt, ...) {
    stop(sprintf("`%s` %s", arg, sprintf(fmt, ...)), call. = FALSE)
  }
  if (is.data.frame(value)) {
    num <- vapply(value, function(col) is.numeric(col) || is.logical(col),
                  logical(1))
    if (!all(num)) {
      bad(paste0("has non-numeric columns: %s. Every column must be numeric; ",
                 "reduce them to numeric summaries first."),
          paste(names(value)[!num], collapse = ", "))
    }
    return(as.matrix(value))
  }
  if (is.list(value)) {
    bad("must be a numeric vector, matrix or data frame, not a list.")
  }
  if (!is.numeric(value) && !is.logical(value)) {
    bad("must be numeric, but it is of type %s.", typeof(value))
  }
  value
}

#' Validate a matrix argument at a public boundary
#'
#' Unlike [as_theta_matrix()], which reshapes whatever it is given, this errors
#' rather than guess. A bare vector is read as a single row and must therefore
#' have exactly `d` entries; a length that does not match used to be recycled
#' into a matrix of the right width, which turned one wrong-length parameter
#' vector into several parameter sets and a plausible-looking answer.
#'
#' @param value The user's value: a numeric vector, matrix or data frame.
#' @param d Required number of columns, or `NULL` to accept any width.
#' @param arg Name of the argument, as it appears in the user's call.
#' @param what Optional phrase describing what a column means, e.g.
#'   `"one parameter per column"`. Shown in parentheses.
#' @return A numeric matrix with `d` columns, column names preserved.
#' @keywords internal
check_matrix <- function(value, d = NULL, arg, what = NULL) {
  bad <- function(fmt, ...) {
    stop(sprintf("`%s` %s", arg, sprintf(fmt, ...)), call. = FALSE)
  }
  detail <- if (is.null(what)) "" else sprintf(" (%s)", what)
  width <- function() {
    if (is.null(d)) "must be a numeric vector or matrix" else
      sprintf("must have %s%s", n_things(d, "column"), detail)
  }

  if (is.null(value)) {
    bad("is missing. It %s.", if (is.null(d)) width() else
      sprintf("must be a numeric vector or matrix with %s%s",
              n_things(d, "column"), detail))
  }
  value <- check_numeric(value, arg)
  if (length(dim(value)) > 2L) {
    bad("must be a vector or a matrix, but it is a %d-dimensional array.",
        length(dim(value)))
  }

  nm <- if (is.null(dim(value))) names(value) else colnames(value)
  if (is.null(dim(value))) {
    if (length(value) == 0L) {
      bad("is empty. It %s.", width())
    }
    if (!is.null(d) && length(value) != d) {
      bad(paste0("%s, but it is a length-%d vector. A bare vector is read as ",
                 "a single row; for several rows pass a matrix with %s."),
          width(), length(value), n_things(d, "column"))
    }
    value <- matrix(value, nrow = 1L)
  } else if (!is.null(d) && ncol(value) != d) {
    # Transposed input is the common way to land here, and the row count says
    # so plainly enough to be worth mentioning.
    hint <- if (nrow(value) == d) " Did you mean to transpose it?" else ""
    bad("%s, but it has %d.%s", width(), ncol(value), hint)
  }

  storage.mode(value) <- "double"
  if (is.null(colnames(value)) && !is.null(nm) && length(nm) == ncol(value)) {
    colnames(value) <- nm
  }
  value
}

#' Validate a count argument
#'
#' One finite whole number, at least `min`. `as.integer()` on its own accepts
#' `2.7` and `TRUE` and turns `NA` into a silent shortfall, so counts that
#' decide how many simulations to run or how many draws to keep go through
#' here instead.
#'
#' @param n The user's value.
#' @param arg Name of the argument.
#' @param min Smallest allowed value.
#' @param why Optional clause explaining the bound, appended to the message.
#' @return `n` as an integer.
#' @keywords internal
check_count <- function(n, arg, min = 1L, why = NULL) {
  ok <- is.numeric(n) && length(n) == 1L && !is.na(n) && is.finite(n) &&
    n == trunc(n) && n >= min
  if (!ok) {
    stop(sprintf("`%s` must be a single whole number of at least %d%s, not %s.",
                 arg, min, if (is.null(why)) "" else paste0(" ", why),
                 describe_value(n)),
         call. = FALSE)
  }
  as.integer(n)
}

#' Validate a probability argument
#'
#' @param p The user's value.
#' @param arg Name of the argument.
#' @param open Require `0 < p < 1` (the default). `FALSE` allows the endpoints.
#' @return `p` as a double.
#' @keywords internal
check_prob <- function(p, arg, open = TRUE) {
  ok <- is.numeric(p) && length(p) == 1L && !is.na(p) && is.finite(p) &&
    (if (isTRUE(open)) p > 0 && p < 1 else p >= 0 && p <= 1)
  if (!ok) {
    stop(sprintf("`%s` must be a single number %s, not %s.", arg,
                 if (isTRUE(open)) "strictly between 0 and 1" else
                   "between 0 and 1 inclusive",
                 describe_value(p)),
         call. = FALSE)
  }
  as.double(p)
}

#' Validate a vector of probabilities
#'
#' [check_prob()] for an argument that is a vector rather than one number, such
#' as the credibility `levels` of [expected_coverage()]. A level outside `(0,
#' 1)` produces one meaningless number in a table of plausible ones:
#' `levels = c(-1, 2)` scores the empty interval and the whole line, and comes
#' back as coverage 0 and 1 with nothing said. The message lists the values, as
#' [check_counts()] does, since which entry is wrong is the thing worth reading.
#'
#' @param p The user's value.
#' @param arg Name of the argument.
#' @param open Require `0 < p < 1` (the default). `FALSE` allows the endpoints.
#' @return `p` as a double vector.
#' @keywords internal
check_probs <- function(p, arg, open = TRUE) {
  ok <- is.numeric(p) && length(p) >= 1L && !anyNA(p) && all(is.finite(p)) &&
    (if (isTRUE(open)) all(p > 0 & p < 1) else all(p >= 0 & p <= 1))
  if (!ok) {
    shown <- if (is.numeric(p) && length(p) > 1L) {
      paste(vapply(p, describe_value, character(1)), collapse = ", ")
    } else {
      describe_value(p)
    }
    stop(sprintf("`%s` must be numbers %s, not %s.", arg,
                 if (isTRUE(open)) "strictly between 0 and 1" else
                   "between 0 and 1 inclusive",
                 shown),
         call. = FALSE)
  }
  as.double(p)
}

#' Validate a column index, which may be given as a name
#'
#' `plot_sbc(sbc_result, param = 99)` used to reach `ranks[, 99]` and report
#' "subscript out of bounds", which names neither the argument nor how many
#' parameters there are. The rank matrix carries `colnames` and every other
#' plotting function labels by name, so a name is accepted here too and
#' resolved against `nms`.
#'
#' @param value The user's value: one index, or one name to match against `nms`.
#' @param arg Name of the argument.
#' @param nms Names to match a character value against, or `NULL` when the
#'   columns are unnamed.
#' @param n Number of columns.
#' @param what Noun for one column, e.g. `"parameter"`.
#' @return An integer index between 1 and `n`.
#' @keywords internal
check_index <- function(value, arg, nms = NULL, n, what = "parameter") {
  bad <- function(fmt, ...) {
    stop(sprintf("`%s` %s", arg, sprintf(fmt, ...)), call. = FALSE)
  }
  if (is.factor(value)) value <- as.character(value)
  if (is.character(value)) {
    if (length(value) != 1L || is.na(value)) {
      bad("must be one %s name or index, not %s.", what, describe_value(value))
    }
    if (is.null(nms)) {
      bad(paste0("is \"%s\", but the %ss are unnamed. Give an index between 1 ",
                 "and %d instead."),
          value, what, n)
    }
    i <- match(value, nms)
    if (is.na(i)) {
      bad("is \"%s\", which is not one of the %s names: %s.",
          value, what, paste(nms, collapse = ", "))
    }
    return(as.integer(i))
  }
  ok <- is.numeric(value) && length(value) == 1L && !is.na(value) &&
    is.finite(value) && value == trunc(value) && value >= 1 && value <= n
  if (!ok) {
    bad("must be one %s index between 1 and %d%s, not %s.", what, n,
        if (is.null(nms)) "" else
          sprintf(", or one of %s", paste(nms, collapse = ", ")),
        describe_value(value))
  }
  as.integer(value)
}

#' Validate a file path argument
#'
#' `saveRDS()` and `writeLines()` take a connection as well as a path, so a
#' value that is neither gets as far as the file system before anything
#' complains. On the reading side `readRDS()` on a path that does not exist
#' raises "cannot open the connection", and the file it could not open is named
#' only in the accompanying warning.
#'
#' @param path The user's value.
#' @param arg Name of the argument, as it appears in the user's call.
#' @param must_exist Require the file to exist, for a path that is read.
#' @return `path`, unchanged.
#' @keywords internal
check_path <- function(path, arg = "path", must_exist = FALSE) {
  ok <- is.character(path) && length(path) == 1L && !is.na(path) &&
    nzchar(path)
  if (!ok) {
    empty <- is.character(path) && length(path) == 1L && !is.na(path)
    stop(sprintf("`%s` must be a single file path, not %s.", arg,
                 if (empty) "an empty string" else describe_value(path)),
         call. = FALSE)
  }
  if (isTRUE(must_exist) && !file.exists(path)) {
    stop(sprintf("`%s` is \"%s\", which does not exist.", arg, path),
         call. = FALSE)
  }
  path
}

#' Validate a vector of counts
#'
#' [check_count()] for an argument that is a vector of counts rather than one,
#' such as `hidden`, where each entry is a layer width and a single bad entry
#' breaks the network the same way a bad scalar would. The message lists the
#' values rather than describing the vector, since which entry is wrong is the
#' thing worth reading.
#'
#' @param n The user's value.
#' @param arg Name of the argument.
#' @param min Smallest allowed value.
#' @param what Optional phrase describing what one entry means, e.g.
#'   `"one hidden-layer width per entry"`. Shown in parentheses.
#' @return `n` as an integer vector.
#' @keywords internal
check_counts <- function(n, arg, min = 1L, what = NULL) {
  ok <- is.numeric(n) && length(n) >= 1L && !anyNA(n) && all(is.finite(n)) &&
    all(n == trunc(n)) && all(n >= min)
  if (!ok) {
    shown <- if (is.numeric(n) && length(n) > 1L) {
      paste(vapply(n, describe_value, character(1)), collapse = ", ")
    } else {
      describe_value(n)
    }
    stop(sprintf("`%s` must be whole numbers of at least %d%s, not %s.",
                 arg, min, if (is.null(what)) "" else sprintf(" (%s)", what),
                 shown),
         call. = FALSE)
  }
  as.integer(n)
}

#' Validate a strictly positive scalar
#'
#' @param x The user's value.
#' @param arg Name of the argument.
#' @param allow_inf Accept `Inf`, for a bound that is disabled by setting it
#'   to infinity (`clip_grad_norm`).
#' @return `x` as a double.
#' @keywords internal
check_positive <- function(x, arg, allow_inf = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && !is.na(x) && x > 0 &&
    (isTRUE(allow_inf) || is.finite(x))
  if (!ok) {
    stop(sprintf("`%s` must be a single positive number%s, not %s.", arg,
                 if (isTRUE(allow_inf)) " or Inf" else "", describe_value(x)),
         call. = FALSE)
  }
  as.double(x)
}

#' Validate a callback argument
#'
#' A function stored now and called later fails at the call site, which can be
#' a training run and several frames away from the argument that was wrong.
#' Checking at the constructor that it is a function, and that it can take the
#' argument it will be given, keeps the complaint next to the mistake.
#'
#' @param f The user's value.
#' @param arg Name of the argument.
#' @param what Optional phrase naming the argument the function receives, e.g.
#'   `"the number of draws"`. Shown in parentheses.
#' @return `f`, invisibly.
#' @keywords internal
check_function <- function(f, arg, what = NULL) {
  detail <- if (is.null(what)) "" else sprintf(" (%s)", what)
  if (!is.function(f)) {
    stop(sprintf("`%s` must be a function of one argument%s, not %s.",
                 arg, detail, describe_value(f)),
         call. = FALSE)
  }
  fmls <- tryCatch(formals(args(f)), error = function(e) NULL)
  if (length(fmls) == 0L) {
    stop(sprintf("`%s` must be a function of one argument%s, but it takes none.",
                 arg, detail),
         call. = FALSE)
  }
  invisible(f)
}

#' Validate a support bound
#'
#' [within_support()] compares `theta` against `lower`/`upper` with `sweep()`,
#' which recycles a bound of the wrong length and warns instead of stopping.
#' The support test that comes back is then wrong, and it decides which
#' posterior draws are rejected as leakage and what `log_prob()` renormalizes
#' by. A length-1 bound is recycled here, once, so everything downstream sees
#' one bound per parameter.
#'
#' @param value The user's value, or `NULL` for an unbounded side.
#' @param arg Name of the argument.
#' @param d Number of parameters.
#' @return `NULL`, or a double vector of length `d`.
#' @keywords internal
check_bound <- function(value, arg, d) {
  if (is.null(value)) return(NULL)
  ok <- is.numeric(value) && (length(value) == 1L || length(value) == d) &&
    !anyNA(value)
  if (!ok) {
    stop(sprintf(paste0("`%s` must be numeric of length %d (one bound per ",
                        "parameter) or length 1, not %s."),
                 arg, d, describe_value(value)),
         call. = FALSE)
  }
  rep_len(as.double(value), d)
}

#' Validate a prior argument, and optionally its dimension
#'
#' @param prior The user's value.
#' @param arg Name of the argument.
#' @param dim Number of parameters the caller expects, or `NULL` to skip that
#'   check.
#' @return The prior, invisibly.
#' @keywords internal
check_prior <- function(prior, arg = "prior", dim = NULL) {
  if (!inherits(prior, "nsbi_prior")) {
    stop(sprintf(paste0("`%s` must be an nsbi_prior object, not %s. Build one ",
                        "with prior_uniform(), prior_normal() or ",
                        "prior_custom(); see ?priors."),
                 arg, if (is.null(prior)) "NULL" else class(prior)[1L]),
         call. = FALSE)
  }
  if (!is.null(dim) && !identical(as.integer(prior$dim), as.integer(dim))) {
    stop(sprintf(paste0("`%s` covers %s but %d %s expected here. The prior ",
                        "must describe exactly the parameters being inferred."),
                 arg, n_things(prior$dim, "parameter"), dim,
                 if (dim == 1L) "is" else "are"),
         call. = FALSE)
  }
  invisible(prior)
}

#' Require every entry to be finite
#'
#' `NA`, `NaN` and `Inf` all reach the estimators as a `chol()` failure or a
#' non-finite validation loss, which blames training for a bad input. Naming
#' the argument and the first offending position at the boundary is cheaper to
#' act on.
#'
#' @param m A numeric vector or matrix.
#' @param arg Name of the argument.
#' @return `m`, invisibly.
#' @keywords internal
check_finite <- function(m, arg) {
  if (!is.numeric(m) && !is.logical(m)) {
    stop(sprintf("`%s` must be numeric, but it is of type %s.", arg, typeof(m)),
         call. = FALSE)
  }
  bad <- !is.finite(m)
  if (!any(bad)) return(invisible(m))
  kinds <- c("NA", "NaN", "Inf")[c(any(is.na(m) & !is.nan(m)), any(is.nan(m)),
                                   any(is.infinite(m)))]
  first <- which(bad)[1L]
  where <- if (is.null(dim(m))) {
    sprintf("position %d", first)
  } else {
    sprintf("row %d, column %d", (first - 1L) %% nrow(m) + 1L,
            (first - 1L) %/% nrow(m) + 1L)
  }
  stop(sprintf("`%s` contains %s (%s), first at %s.",
               arg, n_things(sum(bad), "non-finite value"),
               paste(kinds, collapse = "/"), where),
       call. = FALSE)
}
