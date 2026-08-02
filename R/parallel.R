#' Running the simulator in parallel
#'
#' Every `neuralsbi` function that calls a simulator -- [npe()],
#' [simulate_for_sbi()], [npe_sequential()], [sbc()], [tarp()],
#' [posterior_predictive()] -- runs it through one execution path. That path is
#' sequential by default and parallel as soon as you declare a \pkg{future}
#' plan:
#'
#' ```r
#' library(future)
#' plan(multisession)          # or multicore, cluster, batchtools, ...
#' fit <- npe(prior, simulator, n_simulations = 10000)
#' ```
#'
#' Nothing else changes: no argument to set, no variant of the function to
#' call. With no plan (or `plan(sequential)`) the simulations run in the
#' current process, and `neuralsbi` mentions the two lines above once per
#' session. Silence that hint with
#' `options(neuralsbi.parallel_hint = FALSE)`.
#'
#' @section Random numbers:
#'
#' Each *simulation* gets its own L'Ecuyer-CMRG random-number stream, derived
#' from the session's RNG state at the moment simulation starts. So
#' `set.seed(42)` before a fit (or `npe(..., seed = 42)`) gives the same
#' simulations whatever the plan and whatever the worker count: results depend
#' on the seed alone.
#'
#' Because the simulations are drawn from separate streams, the simulator no
#' longer consumes the caller's RNG state directly; the state advances by one
#' draw per simulation phase.
#'
#' @section Scheduling:
#'
#' Under a multi-process plan the draws are dealt out to workers in batches,
#' because one `future` per simulation would spend longer shipping the
#' simulator to a worker than running it. Batch sizes follow the worker count
#' and are not a setting: they cannot change a result, since every simulation
#' has its own RNG stream. Running sequentially there are no batches at all,
#' just a loop.
#'
#' Each batch ships the simulator, and everything its environment captures, to
#' a worker. A simulator closing over a large object therefore pays for that
#' object once per batch, which is a reason to pass it through `sim_args`
#' instead of capturing it.
#'
#' @section What is not parallelized:
#'
#' Training runs in the calling process. Torch tensors and modules are external
#' pointers that cannot be shipped to a worker, and libtorch already uses
#' multiple threads internally (`torch::torch_set_num_threads()`). The same
#' goes for the posterior-sampling loops in [sbc()] and [tarp()] -- those call
#' the trained network, not the simulator, and are reported with their own
#' progress bar.
#'
#' @name nsbi_parallel
#' @seealso [nsbi_progress] for controlling progress bars.
NULL

#' Number of workers the current future plan provides (1 = sequential)
#' @keywords internal
nsbi_workers <- function() {
  if (!requireNamespace("future", quietly = TRUE)) return(1L)
  n <- tryCatch(as.numeric(future::nbrOfWorkers()), error = function(e) 1)
  if (length(n) != 1L || is.na(n) || n < 1) return(1L)
  # some backends report an unbounded worker count; cap it at something a
  # scheduler can hold in flight
  if (!is.finite(n)) return(128L)
  as.integer(n)
}

#' Tell the user once per session how to go parallel
#' @keywords internal
hint_parallel <- function(n = Inf) {
  if (nsbi_workers() > 1L) return(invisible(FALSE))
  if (!isTRUE(getOption("neuralsbi.parallel_hint", TRUE))) return(invisible(FALSE))
  if (isTRUE(.nsbi_state$parallel_hinted)) return(invisible(FALSE))
  if (n < 100L) return(invisible(FALSE))
  .nsbi_state$parallel_hinted <- TRUE
  message("Running the simulator sequentially. To use all your cores:\n",
          "  library(future)\n",
          "  plan(multisession)\n",
          "Hide this hint with options(neuralsbi.parallel_hint = FALSE).")
  invisible(TRUE)
}

#' Deal `seq_len(n)` out to the current plan's workers
#'
#' Batching is a scheduling detail, not a setting: it exists only because one
#' `future` per simulation would cost more in dispatch than it saves in
#' parallelism, and it cannot change a result because every simulation carries
#' its own RNG stream. A few batches per worker, so uneven simulator run times
#' even out without re-serializing the simulator hundreds of times. Sequential
#' runs get one batch, i.e. a plain loop.
#' @keywords internal
sim_batches <- function(n, workers = nsbi_workers()) {
  n <- as.integer(n)
  if (n <= 0L) return(list())
  if (workers <= 1L) return(list(seq_len(n)))
  size <- ceiling(n / min(n, workers * 4L))
  unname(split(seq_len(n), ceiling(seq_len(n) / size)))
}

#' Independent RNG streams, one per chunk
#'
#' Derived from the caller's current RNG state, so they are reproducible under
#' `set.seed()` yet independent of the backend. Consumes exactly one draw from
#' the caller's stream and leaves its kind untouched.
#' @keywords internal
rng_streams <- function(n) {
  if (n <= 0L) return(list())
  base_seed <- sample.int(.Machine$integer.max, 1L)
  saved <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(assign(".Random.seed", saved, envir = globalenv()), add = TRUE)
  set.seed(base_seed, kind = "L'Ecuyer-CMRG")
  seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  out <- vector("list", n)
  for (i in seq_len(n)) {
    seed <- parallel::nextRNGStream(seed)
    out[[i]] <- seed
  }
  out
}

#' Evaluate `expr` under a given L'Ecuyer-CMRG stream
#'
#' The sequential counterpart of `future(..., seed = stream)`, so sequential
#' and parallel runs draw the same random numbers.
#' @keywords internal
with_rng_stream <- function(stream, expr) {
  saved <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(saved)) {
      suppressWarnings(rm(".Random.seed", envir = globalenv()))
    } else {
      assign(".Random.seed", saved, envir = globalenv())
    }
  }, add = TRUE)
  assign(".Random.seed", stream, envir = globalenv())
  expr
}

#' Apply `fun` to each batch, sequentially or across future workers
#'
#' The parallel branch keeps at most one future per worker in flight and polls
#' for completions, so progress is reported as batches finish rather than in
#' one jump at the end. `fun` is called as `fun(batch, tick)`: running
#' sequentially it ticks the bar itself, once per unit of work; running on a
#' worker it cannot report back mid-batch, so `tick` does nothing and the batch
#' is credited on completion.
#'
#' @param batches List of arguments to `fun`.
#' @param fun Function of a batch and a `tick()` callback.
#' @param p Progress reporter from `nsbi_progressor()`, or `NULL`.
#' @param weights Progress units to credit per batch (defaults to 1 each).
#' @param seeds L'Ecuyer-CMRG seeds, one per batch, handed to the future
#'   backend. `fun` is free to set its own streams inside the batch.
#' @keywords internal
nsbi_batch_apply <- function(batches, fun, p = NULL, weights = NULL,
                             seeds = NULL) {
  n <- length(batches)
  if (n == 0L) return(list())
  weights <- weights %||% rep(1, n)
  seeds <- seeds %||% rng_streams(n)
  workers <- nsbi_workers()

  if (workers <= 1L) {
    tick <- if (is.null(p)) function() NULL else function() p(1)
    out <- vector("list", n)
    for (i in seq_len(n)) {
      out[[i]] <- with_rng_stream(seeds[[i]], fun(batches[[i]], tick))
    }
    return(out)
  }

  noop <- function() NULL
  out <- vector("list", n)
  futures <- vector("list", n)
  pending <- integer(0)   # indices submitted but not yet collected
  nxt <- 1L
  collected <- 0L
  while (collected < n) {
    while (length(pending) < workers && nxt <= n) {
      futures[[nxt]] <- future::futureCall(FUN = fun,
                                           args = list(batches[[nxt]], noop),
                                           seed = seeds[[nxt]])
      pending <- c(pending, nxt)
      nxt <- nxt + 1L
    }
    done <- vapply(pending, function(i) future::resolved(futures[[i]]),
                   logical(1))
    if (!any(done)) {
      Sys.sleep(0.05)
      next
    }
    for (i in pending[done]) {
      out[[i]] <- future::value(futures[[i]])
      futures[i] <- list(NULL)
      collected <- collected + 1L
      if (!is.null(p)) p(weights[i])
    }
    pending <- pending[!done]
  }
  out
}

#' Run a simulator over a parameter matrix
#'
#' The single point through which every simulator call in the package passes:
#' the per-parameter-set contract, the future backend, RNG streams, and the
#' progress bar all live here. See [nsbi_simulator] and [nsbi_parallel].
#'
#' @param simulator The user's simulator, called once per row of `theta`.
#' @param theta Parameter matrix; its `colnames` are the parameter names used
#'   to decide how the simulator receives its arguments.
#' @param sim_args Named list of extra arguments forwarded to every call.
#' @param label Phase name for the progress bar.
#' @param d Expected number of output columns, if known.
#' @return An `n x d` matrix, one row per row of `theta`.
#' @keywords internal
run_simulator <- function(simulator, theta, sim_args = list(),
                          label = "Simulating", d = NULL) {
  theta <- as_theta_matrix(theta)
  call_one <- simulator_caller(simulator, colnames(theta), sim_args)
  n <- nrow(theta)
  if (n == 0L) {
    return(as_theta_matrix(matrix(numeric(0), nrow = 0, ncol = d %||% 1L), d))
  }
  idx <- sim_batches(n)
  hint_parallel(n)
  streams <- rng_streams(n)

  # `finally`, not on.exit(): this block is a promise forced inside
  # progressr::with_progress(), so on.exit() would attach to that frame and
  # close the bar before the work starts.
  with_nsbi_progress({
    p <- nsbi_progressor(steps = n, label = label)
    batches <- lapply(idx, function(i) {
      list(theta = theta[i, , drop = FALSE], index = i, streams = streams[i])
    })
    run_batch <- function(batch, tick) {
      k <- nrow(batch$theta)
      out <- vector("list", k)
      for (j in seq_len(k)) {
        out[[j]] <- as_sim_draw(
          with_rng_stream(
            batch$streams[[j]],
            call_sim_once(call_one, batch$theta[j, ], batch$index[j])
          ),
          batch$index[j]
        )
        tick()
      }
      out
    }
    draws <- tryCatch(
      nsbi_batch_apply(batches, run_batch, p = p, weights = lengths(idx),
                       seeds = lapply(batches, function(b) b$streams[[1L]])),
      finally = p(0, done = TRUE)
    )
    bind_sim_draws(do.call(c, draws), d)
  })
}
