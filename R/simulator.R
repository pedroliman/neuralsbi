#' The simulator contract
#'
#' A simulator in `neuralsbi` is called **once per parameter set** and returns
#' **one simulated observation**. Everything that calls a simulator --
#' [npe()], [simulate_for_sbi()], [npe_sequential()], [sbc()], [tarp()],
#' [posterior_predictive()] -- uses this contract.
#'
#' @section How parameters arrive:
#'
#' Two signatures are accepted, and which one applies is decided once per run
#' from `formals(simulator)`. If every parameter name appears among the
#' simulator's formals, the parameters are passed **by name**, one scalar each:
#'
#' ```r
#' prior <- prior_uniform(low = c(mu = -5, sigma = 0.1),
#'                        high = c(mu =  5, sigma = 3))
#'
#' simulator <- function(mu, sigma) {
#'   y <- rnorm(100, mean = mu, sd = sigma)
#'   c(mean = mean(y), sd = sd(y))
#' }
#' ```
#'
#' Otherwise the whole **named parameter vector** goes to the first argument:
#'
#' ```r
#' simulator <- function(theta) {
#'   y <- rnorm(100, mean = theta["mu"], sd = theta["sigma"])
#'   c(mean = mean(y), sd = sd(y))
#' }
#' ```
#'
#' An unnamed prior always takes the vector form, and so does a prior whose
#' names are not syntactic R names (`"beta[1]"`), since those can never match a
#' formal. A partial match warns and names the parameters that found no formal,
#' because the vector form then hands every parameter to the first argument and
#' the remaining formals fall back to their defaults.
#'
#' @section Everything else the simulator needs:
#'
#' Observed data, a time grid, a population size, a design matrix, a solver
#' tolerance: anything that is not a calibrated parameter goes in `sim_args`, a
#' named list forwarded to every simulator call.
#'
#' ```r
#' simulator <- function(alpha, beta, sigma, x_grid) {
#'   rnorm(length(x_grid), mean = alpha + beta * x_grid, sd = sigma)
#' }
#'
#' fit <- npe(prior, simulator, n_simulations = 10000,
#'            sim_args = list(x_grid = seq(-1, 1, length.out = 50)))
#' ```
#'
#' A list is used rather than `...` because every one of `npe()`'s own
#' arguments would otherwise be exposed to R's partial matching: `x`, `theta`,
#' `n` and `seed` are all natural names for a simulator argument and all of
#' them would be captured silently.
#'
#' @section What the simulator returns:
#'
#' One observation per call: a numeric vector of length `d`, a scalar, or a
#' one-row matrix or data frame. Names on the vector, or `colnames` on the
#' matrix or data frame, become the outcome names used in summaries and plots.
#' Every draw must agree on `d`. Lists, character or factor columns, and
#' anything with more than one row are rejected by name.
#'
#' One sharp edge: arithmetic on a single named parameter keeps the name, so
#' `theta["beta"] / theta["gamma"]` returns a scalar named `beta` and would
#' quietly become an outcome called `beta`. Vectors are unaffected --
#' `theta["alpha"] + theta["beta"] * x` comes back unnamed. Use `unname()`, or
#' name the outcome deliberately with `c(r0 = ...)`.
#'
#' @section Failed simulations:
#'
#' A draw whose output contains `NA`, `NaN` or an infinite value is dropped,
#' together with its parameters, and one warning per run reports the count and
#' the rate. A single non-finite value would otherwise poison the training loss
#' for the whole fit and surface much later as a `NaN` validation loss. The
#' parameters are checked the same way, so a pre-computed pair with a missing
#' `theta` is dropped instead of reaching the density estimator.
#'
#' Dropping conditions on the simulator having succeeded. When failure depends
#' on the parameters -- an ODE that diverges for large `beta`, a model that
#' returns `NA` outside a stability region -- the surviving draws are no longer
#' a sample from the prior, and the fit targets the posterior *given success*
#' rather than the posterior. A high drop rate is a modelling signal, and the
#' honest fix is usually to handle the failure inside the simulator.
#'
#' @name nsbi_simulator
#' @seealso [nsbi_parallel] for how the calls are spread across workers.
NULL

#' How this simulator receives its parameters
#'
#' Decided once per run from `formals()`, never by probing the simulator: a
#' wrong guess would produce a wrong posterior with no error.
#'
#' A partial match warns. When some parameter names appear among the formals
#' and others do not, both signatures are plausible and the vector form wins,
#' which sends the whole parameter vector to the first formal and lets the rest
#' take their defaults. That trains on nonsense without erroring, and one typo
#' in a prior name is enough to cause it. A warning rather than an error,
#' because a vector-signature simulator whose first argument happens to carry a
#' parameter's name is legal.
#' @return `"named"` (one scalar per formal) or `"vector"` (the named parameter
#'   vector as the first argument).
#' @keywords internal
sim_dispatch <- function(simulator, param_names) {
  if (!is.function(simulator)) {
    stop("`simulator` must be a function of one parameter set. See ",
         "?nsbi_simulator.", call. = FALSE)
  }
  if (is.null(param_names) || !all(nzchar(param_names))) return("vector")
  fmls <- names(formals(args(simulator)))
  if (is.null(fmls)) return("vector")
  named <- setdiff(fmls, "...")
  matched <- param_names %in% named
  if (all(matched)) return("named")
  if (any(matched)) {
    warning(sprintf(paste0(
      "Simulator formals match some parameter names but not all: %s %s no ",
      "formal. Passing the whole parameter vector to `%s` instead, which is ",
      "probably not what you want. See ?nsbi_simulator."),
      paste(param_names[!matched], collapse = ", "),
      if (sum(!matched) == 1L) "has" else "have",
      fmls[1L]), call. = FALSE)
  }
  "vector"
}

#' Build the closure that calls the simulator for one parameter set
#'
#' The dispatch decision and the `sim_args` checks happen here, once, so the
#' per-draw loop is a single `do.call()`.
#' @keywords internal
simulator_caller <- function(simulator, param_names, sim_args = list()) {
  mode <- sim_dispatch(simulator, param_names)
  if (length(sim_args)) {
    if (!is.list(sim_args)) {
      stop("`sim_args` must be a named list of extra simulator arguments.",
           call. = FALSE)
    }
    nm <- names(sim_args)
    if (is.null(nm) || !all(nzchar(nm))) {
      stop("Every element of `sim_args` must be named; they are matched to ",
           "the simulator's arguments by name.", call. = FALSE)
    }
    clash <- intersect(nm, param_names %||% character(0))
    if (mode == "named" && length(clash)) {
      stop(sprintf(
        "`sim_args` names clash with parameter names: %s. A name cannot be both a calibrated parameter and a fixed argument.",
        paste(clash, collapse = ", ")), call. = FALSE)
    }
  }
  if (mode == "named") {
    function(theta_i) do.call(simulator, c(as.list(theta_i), sim_args))
  } else {
    function(theta_i) do.call(simulator, c(list(theta_i), sim_args))
  }
}

#' Call the simulator for one parameter set, naming it if it fails
#'
#' [as_sim_draw()] reports everything wrong with a simulator's *return value*
#' and names the simulation it came from. A failure inside the simulator
#' happens before there is a return value, so on its own it arrives as bare as
#' R left it: `argument "ls" is missing, with no default`, with no index, no
#' parameters, and nothing to suggest that the real problem is a dispatch
#' mismatch rather than a missing default. Under a \pkg{future} plan the call
#' crosses a worker boundary first, which is where a useful traceback goes.
#'
#' The original condition is kept as the parent, and its class is carried onto
#' the re-raised error, so a caller catching a condition the simulator signals
#' still catches it.
#'
#' @param call_one The closure from [simulator_caller()].
#' @param theta_i One parameter set.
#' @param i Index of the simulation, for the error message.
#' @return Whatever the simulator returned.
#' @keywords internal
call_sim_once <- function(call_one, theta_i, i = 1L) {
  tryCatch(call_one(theta_i), error = function(e) {
    stop(errorCondition(
      sprintf(paste0("Simulation %d failed: %s\n  parameters: %s\n",
                     "  See ?nsbi_simulator for the two accepted simulator ",
                     "signatures."),
              i, sub("[[:space:]]+$", "", conditionMessage(e)),
              describe_params(theta_i)),
      parent = e,
      class = unique(c("nsbi_sim_error",
                       setdiff(class(e), c("error", "condition"))))))
  })
}

#' Coerce one simulator return value to a named numeric vector
#'
#' Every rejection names the actual problem, because the alternative is a
#' coercion error from three frames down inside `as.matrix()`.
#' @param out What the simulator returned.
#' @param i Index of the simulation, for the error message.
#' @keywords internal
as_sim_draw <- function(out, i = 1L) {
  bad <- function(fmt, ...) {
    stop(sprintf("Simulation %d: %s", i, sprintf(fmt, ...)), call. = FALSE)
  }
  if (is.data.frame(out)) {
    if (nrow(out) != 1L) {
      bad(paste0("the simulator returned a data frame with %d rows. It must ",
                 "return one simulation per call; see ?nsbi_simulator."),
          nrow(out))
    }
    num <- vapply(out, function(col) is.numeric(col) || is.logical(col),
                  logical(1))
    if (!all(num)) {
      bad(paste0("the simulator returned a data frame with non-numeric ",
                 "column(s): %s. Reduce the output to numeric summaries ",
                 "inside the simulator."),
          paste(names(out)[!num], collapse = ", "))
    }
    v <- as.numeric(unlist(out, use.names = FALSE))
    names(v) <- names(out)
    return(v)
  }
  if (!is.null(dim(out))) {
    dm <- dim(out)
    if (length(dm) > 2L) {
      bad(paste0("the simulator returned a %d-dimensional array. It must ",
                 "return a vector, a scalar, or a one-row matrix."), length(dm))
    }
    if (!is.numeric(out) && !is.logical(out)) {
      bad("the simulator returned a %s matrix. The output must be numeric.",
          typeof(out))
    }
    if (dm[1L] != 1L) {
      bad(paste0("the simulator returned a matrix with %d rows. It is called ",
                 "once per parameter set and must return one simulation; see ",
                 "?nsbi_simulator."), dm[1L])
    }
    v <- as.numeric(out)
    names(v) <- colnames(out)
    return(v)
  }
  if (is.list(out)) {
    bad(paste0("the simulator returned a list. It must return a numeric ",
               "vector, a scalar, or a one-row matrix or data frame."))
  }
  if (!is.numeric(out) && !is.logical(out)) {
    bad("the simulator returned %s. The output must be numeric.",
        class(out)[1L])
  }
  if (length(out) == 0L) bad("the simulator returned nothing (length 0).")
  v <- as.numeric(out)
  names(v) <- names(out)
  v
}

#' Stack one-observation draws into an `n x d` matrix
#' @keywords internal
bind_sim_draws <- function(draws, d = NULL) {
  n <- length(draws)
  if (n == 0L) {
    return(as_theta_matrix(matrix(numeric(0), nrow = 0, ncol = d %||% 1L), d))
  }
  lens <- lengths(draws)
  if (any(lens != lens[1L])) {
    j <- which(lens != lens[1L])[1L]
    stop(sprintf(
      paste0("Simulation %d returned %d value(s) but simulation 1 returned ",
             "%d. Every simulation must produce the same number of outputs."),
      j, lens[j], lens[1L]), call. = FALSE)
  }
  x <- matrix(unlist(draws, use.names = FALSE), nrow = n, byrow = TRUE)
  nm <- names(draws[[1L]])
  if (!is.null(nm) && length(nm) == ncol(x) && all(nzchar(nm))) {
    colnames(x) <- nm
  }
  as_theta_matrix(x, d)
}

#' Drop simulations whose parameters or output are not finite
#'
#' A row survives only when every value in both `theta` and `x` is finite, and
#' the two are subset together so the training set stays paired. Checking
#' `theta` matters on the pre-computed path, where a user-supplied `NA` would
#' otherwise reach `chol()` in the estimator and surface as a linear-algebra
#' error. Returns the surviving pair and the number dropped; warns once, and
#' errors when nothing survives. See the "Failed simulations" section of
#' [nsbi_simulator].
#'
#' @param theta Parameter matrix, or `NULL` when the caller has no parameters
#'   to keep aligned.
#' @param x Simulated data matrix.
#' @param what Noun for the warning ("simulations", "SBC trials", ...).
#' @param x_hint What to check when the non-finite values are in `x`. The
#'   default points at the simulator, which is only useful advice when one was
#'   run; the pre-computed `(theta, x)` path passes advice about `x` itself.
#' @keywords internal
drop_failed_sims <- function(theta, x, what = "simulations",
                             x_hint = "Check the simulator on a single prior draw.") {
  n <- nrow(x)
  if (n == 0L) return(list(theta = theta, x = x, n_dropped = 0L, ok = logical(0)))
  bad_x <- rowSums(!is.finite(x)) > 0L
  bad_theta <- if (is.null(theta)) {
    rep(FALSE, n)
  } else {
    rowSums(!is.finite(theta)) > 0L
  }
  ok <- !(bad_x | bad_theta)
  n_dropped <- sum(!ok)
  if (n_dropped > 0L) {
    cause <- if (!any(bad_theta)) {
      "output"
    } else if (!any(bad_x)) {
      "parameters"
    } else {
      "parameters or output"
    }
    if (!any(ok)) {
      hint <- if (any(bad_x)) x_hint else "Check `theta` for NA, NaN or Inf."
      stop(sprintf("All %d %s returned non-finite %s (NA, NaN or Inf). %s",
                   n, what, cause, hint), call. = FALSE)
    }
    warning(sprintf("Dropped %d of %d %s with non-finite %s (%.1f%%).",
                    n_dropped, n, what, cause, 100 * n_dropped / n),
            call. = FALSE)
  }
  list(theta = if (is.null(theta)) NULL else theta[ok, , drop = FALSE],
       x = x[ok, , drop = FALSE], n_dropped = n_dropped, ok = ok)
}
