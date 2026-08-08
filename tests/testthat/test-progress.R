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

test_that("nsbi_progressor() is a no-op under auto mode with no progressr and no interactive session", {
  # progress_backend() falls back to "builtin" here (auto, has_progressr()
  # FALSE), but progress_render() is also FALSE (auto, non-interactive), so
  # nsbi_progressor() must still return the silent no-op reporter rather than
  # drawing the built-in bar.
  opts <- options(neuralsbi.progress = "auto")
  on.exit(options(opts), add = TRUE)
  local_mocked_bindings(has_progressr = function() FALSE)
  expect_equal(progress_backend(), "builtin")
  expect_false(progress_render())

  p <- nsbi_progressor(steps = 5, label = "Testing")
  expect_silent(for (i in 1:5) p(1))
  expect_silent(p(0, done = TRUE))
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

test_that("with_nsbi_progress() wraps expr under its own handler when progressr is active", {
  # with_nsbi_progress() is the wrapper sbc()/tarp()/run_simulator() actually
  # call; the test above exercises nsbi_progressor() directly and so never
  # touches it. Force the "on" option (progress_render() needs it; the
  # non-interactive test session would otherwise read as "auto" == off) and
  # install a silent handler so this doesn't print a bar.
  skip_if_not_installed("progressr")
  opts <- options(neuralsbi.progress = TRUE)
  on.exit(options(opts), add = TRUE)
  expect_equal(progress_backend(), "progressr")
  expect_true(progress_render())

  local_mocked_bindings(default_progressr_handler = progressr::handler_void)
  result <- with_nsbi_progress({
    p <- nsbi_progressor(steps = 3, label = "Testing")
    for (i in 1:3) p(1)
    p(0, done = TRUE)
    42
  })
  expect_equal(result, 42)
})

test_that("with_nsbi_progress() defers to an already-registered global handler", {
  skip_if_not_installed("progressr")
  opts <- options(neuralsbi.progress = TRUE)
  on.exit(options(opts), add = TRUE)
  # Mocking progressr::handlers() rather than actually registering a global
  # handler: real registration is stateful across the whole session and
  # errors if one is already on the calling-handler stack, which makes it
  # fragile to run alongside other tests.
  local_mocked_bindings(handlers = function(...) TRUE, .package = "progressr")

  expect_equal(with_nsbi_progress(1 + 1), 2)
})

test_that("default_progressr_handler() prefers cli, falls back to txtprogressbar", {
  skip_if_not_installed("progressr")
  skip_if_not_installed("cli")
  expect_match(class(default_progressr_handler())[1], "cli")

  local_mocked_bindings(requireNamespace = function(package, ...) !identical(package, "cli"),
                         .package = "base")
  expect_match(class(default_progressr_handler())[1], "txtprogressbar")
})

test_that("builtin_bar() throttles redraws and goes silent once finished", {
  bar <- builtin_bar("Testing")
  out1 <- capture.output(bar(1, 10), type = "message")
  expect_true(any(nzchar(out1)))
  # immediately again: under the 0.1s redraw threshold, so no output
  out2 <- capture.output(bar(2, 10), type = "message")
  expect_length(out2, 0)

  Sys.sleep(0.12)
  out3 <- capture.output(bar(10, 10, done = TRUE), type = "message")
  expect_true(any(nzchar(out3)))
  # finished: further calls are silent even after the throttle window passes
  Sys.sleep(0.12)
  out4 <- capture.output(bar(10, 10, done = TRUE), type = "message")
  expect_length(out4, 0)
})

test_that("nsbi_progressor() grows the total when more progress lands than expected", {
  opts <- options(neuralsbi.progress = "builtin")
  on.exit(options(opts), add = TRUE)
  out <- capture.output({
    p <- nsbi_progressor(steps = 5)
    Sys.sleep(0.12); p(10)   # overshoots the declared 5 steps
    p(0, done = TRUE)
  }, type = "message")
  txt <- paste(out, collapse = "")
  expect_match(txt, "10/10")
})

test_that("durations read as human time", {
  expect_equal(format_duration(0), "0s")
  expect_equal(format_duration(45.4), "45s")
  expect_equal(format_duration(90), "1m30s")
  expect_equal(format_duration(3725), "1h02m")
  expect_equal(format_duration(NA), "?")
  expect_equal(format_duration(Inf), "?")
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
