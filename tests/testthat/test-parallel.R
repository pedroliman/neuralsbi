test_that("batching is scheduling only: sequential runs are one loop", {
  # no plan declared, so no batching at all
  expect_length(sim_batches(1000), 1L)
  expect_equal(unlist(sim_batches(1000)), seq_len(1000))
  expect_length(sim_batches(0), 0L)

  # under a plan, a few batches per worker, and every draw appears once
  idx <- sim_batches(1000, workers = 4L)
  expect_length(idx, 16L)
  expect_equal(unlist(idx), seq_len(1000))

  # fewer draws than batch slots: one draw each, never an empty batch
  idx <- sim_batches(3, workers = 8L)
  expect_length(idx, 3L)
  expect_true(all(lengths(idx) == 1L))
})

test_that("simulations are reproducible under a seed", {
  prior <- prior_uniform(c(-1, -1), c(1, 1))
  sim <- function(theta) theta + stats::rnorm(length(theta), sd = 0.2)
  a <- simulate_for_sbi(sim, prior, 200, seed = 42)
  b <- simulate_for_sbi(sim, prior, 200, seed = 42)
  expect_identical(a, b)

  # each draw carries its own stream, so the same theta always gives the same
  # x whatever else is in the run: re-simulating one draw's parameters on its
  # own stream reproduces that row exactly
  theta1 <- a$theta[1, , drop = FALSE]
  set.seed(42)
  invisible(sample_prior(prior, 200))       # line the streams up as the run did
  expect_identical(run_simulator(sim, theta1)[1, ], a$x[1, ])
})

test_that("each simulation phase draws fresh randomness", {
  prior <- prior_uniform(-1, 1)
  sim <- function(theta) theta + stats::rnorm(1)
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

test_that("accepted output shapes all give one row per simulation", {
  prior <- prior_uniform(-1, 1)

  # a scalar return is one observation, whatever the chunking
  x <- simulate_for_sbi(function(theta) theta^2, prior, 100)$x
  expect_equal(dim(x), c(100L, 1L))

  # a one-row data frame names the outcomes
  x <- simulate_for_sbi(
    function(theta) data.frame(a = theta, b = theta^2), prior, 60)$x
  expect_equal(dim(x), c(60L, 2L))
  expect_equal(colnames(x), c("a", "b"))

  # so does a one-row matrix, and a named vector
  x <- simulate_for_sbi(
    function(theta) matrix(c(theta, theta^2), nrow = 1,
                           dimnames = list(NULL, c("a", "b"))), prior, 20)$x
  expect_equal(colnames(x), c("a", "b"))
  x <- simulate_for_sbi(function(theta) c(a = 1, b = 2), prior, 20)$x
  expect_equal(colnames(x), c("a", "b"))

  # a length-3 vector is one observation of three outcomes, never three draws
  x <- run_simulator(function(theta) c(1, 2, 3), matrix(0, nrow = 1, ncol = 1))
  expect_equal(dim(x), c(1L, 3L))
})

test_that("a simulator written for one parameter set is not silently accepted", {
  prior <- prior_uniform(c(-1, -1), c(1, 1))

  # the 0.4.0 regression: with 2-row chunks this used to pass as 100 draws of
  # one outcome, because the return length happened to match the chunk size
  x <- simulate_for_sbi(function(theta) c(theta[1] * 2, theta[2] * 2),
                        prior, 100)$x
  expect_equal(dim(x), c(100L, 2L))

  # output that changes length between draws is an error naming both draws
  expect_error(
    run_simulator(function(theta) if (theta[1] > 0) c(1, 1) else 1,
                  matrix(c(-1, 1), ncol = 1)),
    "returned 2 value\\(s\\) but simulation 1 returned 1"
  )
})

test_that("output that is not one numeric observation is rejected by name", {
  th <- matrix(0, nrow = 2, ncol = 1)
  expect_error(run_simulator(function(theta) list(1, 2), th),
               "returned a list")
  expect_error(run_simulator(function(theta) "a", th),
               "returned character")
  expect_error(run_simulator(function(theta) matrix(0, nrow = 3, ncol = 2), th),
               "matrix with 3 rows")
  expect_error(run_simulator(function(theta) data.frame(a = 1, b = "x"), th),
               "non-numeric column\\(s\\): b")
  expect_error(run_simulator(function(theta) array(0, c(2, 2, 2)), th),
               "3-dimensional array")
  expect_error(run_simulator(function(theta) numeric(0), th),
               "returned nothing")
})

test_that("non-finite simulations are dropped in pairs, with a warning", {
  prior <- prior_uniform(c(-1, -1), c(1, 1))
  # every draw with a negative first parameter fails
  sim <- function(theta) if (theta[1] < 0) c(NA, NA) else unname(theta)

  expect_warning(sims <- simulate_for_sbi(sim, prior, 200, seed = 3),
                 "non-finite output")
  expect_equal(nrow(sims$theta), nrow(sims$x))
  expect_equal(sims$n_dropped, 200L - nrow(sims$x))
  expect_true(all(is.finite(sims$x)))
  expect_true(all(sims$theta[, 1] >= 0))

  # nothing left is an error, not a warning
  expect_error(
    simulate_for_sbi(function(theta) c(NA_real_, NA_real_), prior, 20),
    "All 20 simulations returned non-finite output"
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
  sim <- function(theta) theta + stats::rnorm(length(theta), sd = 0.2)
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

test_that("a failed trial lowers the effective n_sbc and is reported", {
  set.seed(21)
  prior <- prior_normal(mean = c(0, 0), sd = 1)
  sim <- function(theta) unname(theta) + rnorm(2, sd = 0.4)
  fit <- npe(prior, sim, n_simulations = 1000,
             density_estimator = "linear_gaussian")

  # fail on trials whose first parameter is positive
  flaky <- function(theta) if (theta[1] > 0) c(NA, NA) else unname(theta)
  expect_warning(res <- sbc(fit, flaky, n_sbc = 60, n_posterior_samples = 50,
                            seed = 4),
                 "SBC trials with non-finite output")
  expect_equal(nrow(res$ranks), res$n_sbc)
  expect_lt(res$n_sbc, 60L)
  expect_gt(res$n_dropped, 0L)
  expect_output(print(res), "dropped for non-finite simulator output")
})
