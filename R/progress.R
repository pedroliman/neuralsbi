#' Progress reporting
#'
#' Long-running steps -- running the simulator, training a neural density
#' estimator -- report progress through one mechanism, so a fit that simulates
#' and then trains shows the same style of bar for both phases. Each bar
#' reports the elapsed time and an ETA extrapolated from the work done so far.
#'
#' If the \pkg{progressr} package is installed, `neuralsbi` emits standard
#' progressr updates, so `progressr::handlers()` and
#' `progressr::with_progress()` control reporting exactly as they do for any
#' other progressr-aware package:
#'
#' ```r
#' library(progressr)
#' handlers(global = TRUE)
#' handlers("cli")            # or "progress", "txtprogressbar", "rstudio"
#' fit <- npe(prior, simulator, n_simulations = 10000)
#' ```
#'
#' Without progressr installed (or when it is installed but no handler is
#' registered) `neuralsbi` draws its own bar, which needs no extra packages.
#'
#' Reporting is controlled by `options(neuralsbi.progress = )`:
#'
#' * `"auto"` (default) -- report in interactive sessions, stay quiet in
#'   scripts, `R CMD check`, and knitr.
#' * `TRUE` -- always report, interactive or not.
#' * `FALSE` -- never report.
#' * `"builtin"` -- always report with the built-in bar, ignoring progressr.
#'
#' @section Training progress:
#'
#' While a neural estimator trains, one progress step is one epoch. The total
#' is not known in advance because training stops early, so the bar targets the
#' epoch at which training would stop if the validation loss never improved
#' again (`best epoch + patience`). The target moves out whenever the
#' validation loss improves, so the ETA is a running estimate, not a promise.
#'
#' @name nsbi_progress
#' @seealso [nsbi_parallel] for running the simulator across cores.
NULL

# Package-level session state (progress/parallel hints).
.nsbi_state <- new.env(parent = emptyenv())

#' Resolve the `neuralsbi.progress` option to one of
#' "auto", "on", "off", "builtin", "progressr"
#' @keywords internal
progress_option <- function() {
  opt <- getOption("neuralsbi.progress", "auto")
  if (is.logical(opt) && length(opt) == 1L && !is.na(opt)) {
    opt <- if (opt) "on" else "off"
  }
  opt <- as.character(opt)[1]
  if (identical(opt, "none")) opt <- "off"
  if (!opt %in% c("auto", "on", "off", "builtin", "progressr")) opt <- "auto"
  opt
}

#' @keywords internal
has_progressr <- function() requireNamespace("progressr", quietly = TRUE)

#' Which reporter carries the updates: "progressr", "builtin", or "none"
#' @keywords internal
progress_backend <- function() {
  opt <- progress_option()
  if (opt == "off") return("none")
  if (opt == "builtin") return("builtin")
  if (has_progressr()) "progressr" else "builtin"
}

#' Should *we* draw a bar?
#'
#' Distinct from `progress_backend()`: progressr users who registered their own
#' handlers always get updates, whatever this returns. This governs only
#' whether `neuralsbi` installs a handler itself (or draws its own bar).
#' @keywords internal
progress_render <- function() {
  switch(progress_option(),
    off = FALSE,
    auto = interactive() && !isTRUE(getOption("knitr.in.progress")),
    TRUE
  )
}

#' Report progress for one phase of work
#'
#' Wraps a phase (simulation, training) so the built-in bar or a progressr
#' handler is active for its duration. `expr` is passed on unevaluated so
#' progressr can catch the progression conditions it signals.
#' @keywords internal
with_nsbi_progress <- function(expr) {
  if (progress_backend() != "progressr" || !progress_render()) {
    return(expr)
  }
  # A global handler means the user is in charge; don't nest our own.
  if (isTRUE(tryCatch(progressr::handlers(global = NA), error = function(e) FALSE))) {
    return(expr)
  }
  progressr::with_progress(expr, handlers = default_progressr_handler(),
                           enable = TRUE)
}

#' Pick a progressr handler that can show an ETA
#'
#' `handler_txtprogressbar` is the progressr default and reports no ETA, so
#' prefer the cli handler when cli is available. A user who has registered
#' handlers of their own never gets here.
#' @keywords internal
default_progressr_handler <- function() {
  if (requireNamespace("cli", quietly = TRUE)) {
    return(progressr::handler_cli(enable = TRUE))
  }
  progressr::handler_txtprogressbar(enable = TRUE)
}

#' Create a progress reporter for one phase
#'
#' Returns `function(amount = 1, total = NULL, done = FALSE)`. `total` revises
#' the denominator mid-flight, which is what training needs: the epoch a run
#' will stop at is only known once it stops. Under progressr -- whose step
#' count is fixed at creation -- the revised fraction is mapped back onto the
#' original `steps`, so both reporters show the same picture.
#'
#' @param steps Expected number of units of work.
#' @param label Short phase name shown next to the bar.
#' @keywords internal
nsbi_progressor <- function(steps, label = NULL) {
  steps <- as.numeric(steps)[1]
  backend <- progress_backend()
  if (backend == "none" || !is.finite(steps) || steps <= 0) {
    return(function(amount = 1, total = NULL, done = FALSE) invisible(NULL))
  }

  state <- new.env(parent = emptyenv())
  state$pos <- 0
  state$total <- steps
  bump <- function(amount, total) {
    if (!is.null(total) && is.finite(total) && total > 0) state$total <- total
    state$pos <- state$pos + amount
    if (state$pos > state$total) state$total <- state$pos
  }

  if (backend == "progressr") {
    # on_exit = FALSE: progressr would otherwise finish the progressor when
    # *this* frame exits, which is immediately. The caller closes it instead,
    # and hitting `steps` auto-finishes it.
    p <- progressr::progressor(steps = steps, on_exit = FALSE)
    state$emitted <- 0
    return(function(amount = 1, total = NULL, done = FALSE) {
      bump(amount, total)
      target <- if (done) steps else floor(state$pos / state$total * steps)
      if (target > state$emitted) {
        p(amount = target - state$emitted, message = label)
        state$emitted <- target
      }
      invisible(NULL)
    })
  }

  if (!progress_render()) {
    return(function(amount = 1, total = NULL, done = FALSE) invisible(NULL))
  }
  bar <- builtin_bar(label)
  function(amount = 1, total = NULL, done = FALSE) {
    bump(amount, total)
    bar(state$pos, state$total, done)
  }
}

#' Minimal dependency-free progress bar with an ETA
#'
#' Draws to `stderr()`, redrawing in place. Updates are throttled so a fast
#' inner loop does not spend its time formatting text.
#' @keywords internal
builtin_bar <- function(label = NULL, width = 22L) {
  start <- Sys.time()
  last_draw <- -Inf
  finished <- FALSE
  label <- label %||% ""

  function(pos, total, done = FALSE) {
    if (finished) return(invisible(NULL))
    now <- as.numeric(Sys.time())
    if (!done && now - last_draw < 0.1) return(invisible(NULL))
    last_draw <<- now
    frac <- if (total > 0) min(1, pos / total) else 0
    if (done) frac <- 1
    filled <- round(frac * width)
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    eta <- if (frac > 0 && !done) elapsed * (1 - frac) / frac else 0
    cat(sprintf("\r%-10s [%s%s] %3.0f%% | %s/%s | %s elapsed | ETA %s   ",
                label,
                strrep("=", filled), strrep(" ", width - filled),
                100 * frac, format_count(pos), format_count(total),
                format_duration(elapsed), format_duration(eta)),
        file = stderr())
    if (done) {
      cat("\n", file = stderr())
      finished <<- TRUE
    }
    invisible(NULL)
  }
}

#' @keywords internal
format_count <- function(x) format(round(x), trim = TRUE, scientific = FALSE)

#' @keywords internal
format_duration <- function(secs) {
  if (!is.finite(secs)) return("?")
  secs <- max(0, round(secs))
  if (secs < 60) return(sprintf("%ds", secs))
  if (secs < 3600) return(sprintf("%dm%02ds", secs %/% 60, secs %% 60))
  sprintf("%dh%02dm", secs %/% 3600, (secs %% 3600) %/% 60)
}
