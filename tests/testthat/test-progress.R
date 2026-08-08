test_that("the progress option resolves as documented", {
  opts <- options(neuralsbi.progress = NULL)
  on.exit(options(opts), add = TRUE)

  expect_equal(progress_option(), "auto")
  options(neuralsbi.progress = TRUE);     expect_equal(progress_option(), "on")
  options(neuralsbi.progress = FALSE);    expect_equal(progress_option(), "off")
  options(neuralsbi.progress = "none");   expect_equal(progress_option(), "off")
  options(neuralsbi.progress = "builtin");expect_equal(progress_option(), "builtin")
  options(neuralsbi.progress = "nonsense"); expect_equal(progress_option(), "auto")

  options(neuralsbi.progress = FALSE)
  expect_equal(progress_backend(), "none")
  expect_false(progress_render())
  options(neuralsbi.progress = "builtin")
  expect_equal(progress_backend(), "builtin")
  expect_true(progress_render())
})

test_that("reporting off means no output and no bookkeeping", {
  opts <- options(neuralsbi.progress = FALSE)
  on.exit(options(opts), add = TRUE)
  p <- nsbi_progressor(steps = 10, label = "Testing")
  expect_silent(for (i in 1:10) p(1))
  expect_silent(p(0, done = TRUE))
})

test_that("the built-in bar reports counts, elapsed time and an ETA", {
  opts <- options(neuralsbi.progress = "builtin")
  on.exit(options(opts), add = TRUE)

  out <- capture.output({
    p <- nsbi_progressor(steps = 4, label = "Testing")
    for (i in 1:4) {
      Sys.sleep(0.12)
      p(1)
    }
    p(0, done = TRUE)
  }, type = "message")
  txt <- paste(out, collapse = "")

  expect_match(txt, "Testing")
  expect_match(txt, "1/4")
  expect_match(txt, "4/4")
  expect_match(txt, "100%")
  expect_match(txt, "ETA")
})

test_that("a revised total moves the denominator, never the position back", {
  opts <- options(neuralsbi.progress = "builtin")
  on.exit(options(opts), add = TRUE)

  out <- capture.output({
    p <- nsbi_progressor(steps = 100, label = "Training")
    Sys.sleep(0.12); p(1, total = 40)
    Sys.sleep(0.12); p(1, total = 60)
    p(0, done = TRUE)
  }, type = "message")
  txt <- paste(out, collapse = "")
  expect_match(txt, "1/40")
  expect_match(txt, "2/60")
})

test_that("progressr carries the updates when it is installed", {
  skip_if_not_installed("progressr")
  opts <- options(neuralsbi.progress = "auto")
  on.exit(options(opts), add = TRUE)
  expect_equal(progress_backend(), "progressr")

  # Count the progression conditions a phase signals. The handler goes inside
  # with_progress(), which muffles anything it has already handled.
  seen <- 0L
  progressr::with_progress(
    withCallingHandlers({
      p <- nsbi_progressor(steps = 10, label = "Testing")
      for (i in 1:10) p(1)
      p(0, done = TRUE)
    }, progression = function(c) seen <<- seen + 1L),
    handlers = progressr::handler_void(),
    enable = TRUE   # progressr reports nothing in a non-interactive session
  )
  expect_gte(seen, 10L)
})

test_that("durations read as human time", {
  expect_equal(format_duration(0), "0s")
  expect_equal(format_duration(45.4), "45s")
  expect_equal(format_duration(90), "1m30s")
  expect_equal(format_duration(3725), "1h02m")
  expect_equal(format_duration(NA), "?")
  expect_equal(format_duration(Inf), "?")
})

test_that("bar_rate prefers the recent window over the lifetime average", {
  # A slow first step (setup cost) followed by a steady, faster pace: the
  # lifetime average (elapsed/pos) is dragged down by the slow start, but a
  # window that excludes it should read the true, faster current rate.
  times <- c(6, 7, 8, 9)      # 1s/step over the trimmed window
  positions <- c(3, 4, 5, 6)
  lifetime_rate <- 6 / 9      # what elapsed/pos would have said
  window_rate <- bar_rate(times, positions, elapsed = 9, pos = 6)
  expect_equal(window_rate, 1)
  expect_gt(window_rate, lifetime_rate)
})

test_that("bar_rate falls back to the lifetime average with a thin window", {
  expect_equal(bar_rate(numeric(0), numeric(0), elapsed = 10, pos = 5), 0.5)
  expect_equal(bar_rate(5, 2, elapsed = 5, pos = 2), 0.4)
  # two samples with no elapsed time between them (throttled to the same
  # tick): can't divide by zero, fall back
  expect_equal(bar_rate(c(3, 3), c(1, 2), elapsed = 3, pos = 2), 2 / 3)
  # nothing to go on at all
  expect_true(is.na(bar_rate(numeric(0), numeric(0), elapsed = 0, pos = 0)))
})

test_that("bar_eta converts a rate into seconds remaining", {
  expect_equal(bar_eta(rate = 2, pos = 4, total = 10), 3)
  expect_equal(bar_eta(rate = 1, pos = 10, total = 10), 0)
  expect_equal(bar_eta(rate = 1, pos = 12, total = 10), 0) # never negative
  expect_true(is.na(bar_eta(rate = NA_real_, pos = 4, total = 10)))
  expect_true(is.na(bar_eta(rate = 0, pos = 4, total = 10)))
})

test_that("the training bar projects the early-stopping epoch", {
  # nothing finished yet: the projection is this restart's stopping epoch
  expect_equal(
    train_progress_total(integer(0), best_epoch = 10, patience = 20,
                         max_epochs = 2000, restart = 1, n_restarts = 1),
    30
  )
  # never better than max_epochs
  expect_equal(
    train_progress_total(integer(0), best_epoch = 1990, patience = 20,
                         max_epochs = 2000, restart = 1, n_restarts = 1),
    2000
  )
  # restarts still to come are budgeted at the mean of those already run
  expect_equal(
    train_progress_total(c(100L, 200L), best_epoch = 10, patience = 20,
                         max_epochs = 2000, restart = 3, n_restarts = 4),
    100 + 200 + 30 + 150
  )
})
