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
#' @section Chunking:
#'
#' The parameter matrix is split into row chunks and one chunk at a time is
#' handed to the simulator, so a vectorized simulator keeps most of its
#' vectorization while progress can still be reported and work spread across
#' workers. The number of chunks depends only on the number of simulations
#' (`ceiling(n / chunk_size)`, with `chunk_size` chosen to give about 64
#' chunks), *not* on the number of workers -- which is what makes results
#' reproducible across plans. Set it yourself with the `chunk_size` argument,
#' or globally with `options(neuralsbi.chunks = )`. Larger chunks mean less
#' overhead per call and a coarser progress bar; smaller chunks balance
#' uneven simulator run times better. Under a multi-process plan each chunk
#' ships the simulator -- and everything its environment captures -- to a
#' worker, so a simulator closing over a large object is a reason to raise
#' `chunk_size`.
#'
#' @section Random numbers:
#'
#' Each chunk gets its own L'Ecuyer-CMRG random-number stream, derived from the
#' session's RNG state at the moment simulation starts. So `set.seed(42)`
#' before a fit (or `npe(..., seed = 42)`) gives the same simulations whether
#' you run sequentially or on 32 workers, and whatever the worker count. What
#' does change results is the chunking: fixing `chunk_size` fixes the draws.
#'
#' Because the simulations are drawn from separate streams, the simulator no
#' longer consumes the caller's RNG state directly; the state advances by one
#' draw per simulation phase.
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

#' Split `seq_len(n)` into row chunks
#'
#' The chunk count depends on `n` alone, so the same seed gives the same
#' simulations under any plan. See the "Chunking" section of [nsbi_parallel].
#' @keywords internal
chunk_index <- function(n, chunk_size = NULL) {
  n <- as.integer(n)
  if (n <= 0L) return(list())
  if (is.null(chunk_size)) {
    target <- as.integer(getOption("neuralsbi.chunks", 64L))
    if (is.na(target) || target < 1L) target <- 1L
    chunk_size <- ceiling(n / target)
  }
  chunk_size <- max(1L, as.integer(chunk_size))
  unname(split(seq_len(n), ceiling(seq_len(n) / chunk_size)))
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

#' Apply `fun` to each chunk, sequentially or across future workers
#'
#' The parallel branch keeps at most one future per worker in flight and polls
#' for completions, so progress is reported as chunks finish rather than in one
#' jump at the end.
#'
#' @param chunks List of arguments to `fun`.
#' @param fun Function of one argument.
#' @param p Progress reporter from `nsbi_progressor()`, or `NULL`.
#' @param weights Progress units to credit per chunk (defaults to 1 each).
#' @keywords internal
nsbi_chunk_apply <- function(chunks, fun, p = NULL, weights = NULL) {
  n <- length(chunks)
  if (n == 0L) return(list())
  weights <- weights %||% rep(1, n)
  streams <- rng_streams(n)
  workers <- nsbi_workers()
  tick <- function(i) if (!is.null(p)) p(weights[i])

  if (workers <= 1L) {
    out <- vector("list", n)
    for (i in seq_len(n)) {
      out[[i]] <- with_rng_stream(streams[[i]], fun(chunks[[i]]))
      tick(i)
    }
    return(out)
  }

  out <- vector("list", n)
  futures <- vector("list", n)
  pending <- integer(0)   # indices submitted but not yet collected
  nxt <- 1L
  collected <- 0L
  while (collected < n) {
    while (length(pending) < workers && nxt <= n) {
      futures[[nxt]] <- future::futureCall(FUN = fun, args = list(chunks[[nxt]]),
                                           seed = streams[[nxt]])
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
      tick(i)
    }
    pending <- pending[!done]
  }
  out
}

#' Run a simulator over a parameter matrix
#'
#' The single point through which every simulator call in the package passes:
#' chunking, the future backend, RNG streams, and the progress bar all live
#' here. See [nsbi_parallel].
#'
#' @param simulator The user's simulator.
#' @param theta Parameter matrix.
#' @param chunk_size Rows per simulator call; `NULL` picks a default.
#' @param label Phase name for the progress bar.
#' @param d Expected number of output columns, if known.
#' @keywords internal
run_simulator <- function(simulator, theta, chunk_size = NULL,
                          label = "Simulating", d = NULL) {
  if (!is.function(simulator)) {
    stop("`simulator` must be a function of a parameter matrix.", call. = FALSE)
  }
  theta <- as_theta_matrix(theta)
  n <- nrow(theta)
  if (n == 0L) return(as_theta_matrix(matrix(numeric(0), nrow = 0, ncol = d %||% 1L)))
  idx <- chunk_index(n, chunk_size)
  hint_parallel(n)

  # `finally`, not on.exit(): this block is a promise forced inside
  # progressr::with_progress(), so on.exit() would attach to that frame and
  # close the bar before the work starts.
  with_nsbi_progress({
    p <- nsbi_progressor(steps = n, label = label)
    chunks <- lapply(idx, function(i) theta[i, , drop = FALSE])
    out <- tryCatch(
      nsbi_chunk_apply(chunks, simulator, p = p, weights = lengths(idx)),
      finally = p(0, done = TRUE)
    )
    bind_sim_output(out, lengths(idx), n, d)
  })
}

#' Stack per-chunk simulator output into one matrix
#'
#' A chunk of one row that comes back as a plain vector is a single
#' observation, not a column of them -- the one place where chunking could
#' silently transpose a result.
#' @keywords internal
bind_sim_output <- function(out, sizes, n, d = NULL) {
  mats <- Map(function(v, k) {
    if (is.data.frame(v)) v <- as.matrix(v)
    if (is.null(dim(v))) {
      v <- if (k == 1L) matrix(v, nrow = 1L) else matrix(v, ncol = 1L)
    }
    storage.mode(v) <- "double"
    v
  }, out, sizes)
  ncols <- vapply(mats, ncol, integer(1))
  if (length(unique(ncols)) > 1L) {
    stop("Simulator returned a different number of columns for different ",
         "parameter chunks.", call. = FALSE)
  }
  x <- do.call(rbind, mats)
  if (nrow(x) != n) {
    stop("Simulator must return one row of output per row of theta.",
         call. = FALSE)
  }
  as_theta_matrix(x, d)
}
