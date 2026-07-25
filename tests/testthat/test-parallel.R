test_that("chunking depends on n alone, not on the backend", {
  # ~64 chunks by default, so the same n always splits the same way
  expect_length(chunk_index(1000), 63L)
  expect_equal(lengths(chunk_index(1000))[[1]], 16L)
  expect_length(chunk_index(10), 10L)
  expect_length(chunk_index(1), 1L)
  expect_length(chunk_index(0), 0L)

  idx <- chunk_index(100, chunk_size = 25)
  expect_length(idx, 4L)
  expect_equal(unlist(idx), seq_len(100))

  # a chunk_size larger than n gives a single call, i.e. the old behaviour
  expect_length(chunk_index(100, chunk_size = 1000), 1L)

  withr_option <- options(neuralsbi.chunks = 4L)
  on.exit(options(withr_option), add = TRUE)
  expect_length(chunk_index(100), 4L)
})

test_that("simulations are reproducible under a seed", {
  prior <- prior_uniform(c(-1, -1), c(1, 1))
  sim <- function(theta) theta + matrix(stats::rnorm(length(theta), sd = 0.2),
                                        nrow = nrow(theta))
  a <- simulate_for_sbi(sim, prior, 200, seed = 42)
  b <- simulate_for_sbi(sim, prior, 200, seed = 42)
  expect_identical(a, b)

  # prior draws do not depend on the chunking; only the simulator's own noise
  c1 <- simulate_for_sbi(sim, prior, 200, seed = 42, chunk_size = 200)
  expect_identical(a$theta, c1$theta)
})

test_that("each simulation phase draws fresh randomness", {
  prior <- prior_uniform(-1, 1)
  sim <- function(theta) theta + stats::rnorm(nrow(theta))
  set.seed(1)
  x1 <- simulate_for_sbi(sim, prior, 100)$x
  x2 <- simulate_for_sbi(sim, prior, 100)$x
  expect_false(identical(x1, x2))
})

test_that("the caller's RNG kind survives a simulation", {
  prior <- prior_uniform(-1, 1)
  set.seed(1)
  kind_before <- RNGkind()
  invisible(simulate_for_sbi(function(theta) theta, prior, 50))
  expect_identical(RNGkind(), kind_before)
})

test_that("simulator output shapes survive chunking", {
  prior <- prior_uniform(-1, 1)

  # a plain vector return is one column, not one row
  x <- simulate_for_sbi(function(theta) as.numeric(theta[, 1]^2), prior, 100)$x
  expect_equal(dim(x), c(100L, 1L))

  # a data frame keeps its column names
  x <- simulate_for_sbi(
    function(theta) data.frame(a = theta[, 1], b = theta[, 1]^2), prior, 60)$x
  expect_equal(dim(x), c(60L, 2L))
  expect_equal(colnames(x), c("a", "b"))

  # a one-row chunk returning a vector is a single observation
  x <- run_simulator(function(theta) c(1, 2, 3), matrix(0, nrow = 1, ncol = 1))
  expect_equal(dim(x), c(1L, 3L))
})

test_that("a simulator with the wrong output length is rejected", {
  prior <- prior_uniform(-1, 1)
  expect_error(
    simulate_for_sbi(function(theta) rep(0, nrow(theta) + 1L), prior, 100,
                     chunk_size = 100),
    "one row of output per row of theta"
  )
  expect_error(
    run_simulator(function(theta) {
      if (theta[1, 1] > 0) cbind(1, 1) else 1
    }, matrix(c(-1, 1), ncol = 1), chunk_size = 1),
    "different number of columns"
  )
})

test_that("the parallel hint fires once and can be silenced", {
  old <- .nsbi_state$parallel_hinted
  on.exit(.nsbi_state$parallel_hinted <- old, add = TRUE)

  .nsbi_state$parallel_hinted <- NULL
  expect_message(hint_parallel(1000), "plan\\(multisession\\)")
  expect_silent(hint_parallel(1000))

  .nsbi_state$parallel_hinted <- NULL
  expect_silent(hint_parallel(10))  # too small to be worth parallelizing

  opts <- options(neuralsbi.parallel_hint = FALSE)
  on.exit(options(opts), add = TRUE)
  .nsbi_state$parallel_hinted <- NULL
  expect_silent(hint_parallel(1000))
})

test_that("a future plan gives the same simulations as running sequentially", {
  skip_on_cran()
  skip_if_not_installed("future")

  prior <- prior_uniform(c(-1, -1), c(1, 1))
  sim <- function(theta) theta + matrix(stats::rnorm(length(theta), sd = 0.2),
                                        nrow = nrow(theta))
  seq_run <- simulate_for_sbi(sim, prior, 200, seed = 7)

  old_plan <- future::plan(future::multisession, workers = 2L)
  on.exit(future::plan(old_plan), add = TRUE)
  expect_gt(nsbi_workers(), 1L)
  par_run <- simulate_for_sbi(sim, prior, 200, seed = 7)

  expect_identical(seq_run$theta, par_run$theta)
  expect_identical(seq_run$x, par_run$x)
})

test_that("errors inside a worker reach the caller", {
  skip_on_cran()
  skip_if_not_installed("future")

  prior <- prior_uniform(-1, 1)
  old_plan <- future::plan(future::multisession, workers = 2L)
  on.exit(future::plan(old_plan), add = TRUE)
  expect_error(
    simulate_for_sbi(function(theta) stop("simulator blew up"), prior, 20),
    "simulator blew up"
  )
})
