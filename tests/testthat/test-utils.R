# The Suggests guards. torch and the plotting stack are optional, so the code
# that needs them has to say what to install instead of failing somewhere
# inside a package that is not there. These messages are what a user without
# libtorch sees first, so they are worth pinning.

test_that("require_torch() names what to install and the torch-free fallback", {
  skip_if(torch_available(), "torch works here, so the failing branch is unreachable")
  if (requireNamespace("torch", quietly = TRUE)) {
    expect_error(require_torch(), "libtorch is not\\. Run torch::install_torch")
  } else {
    expect_error(require_torch(), "needs the 'torch' package")
    expect_error(require_torch(), "install\\.packages\\('torch'\\)")
    expect_error(require_torch(), "density_estimator = 'linear_gaussian'")
  }
})

test_that("require_torch() passes when torch is available", {
  skip_if_no_torch()
  expect_true(require_torch())
})

# The two require_torch() error branches above only run in an environment
# that actually lacks torch or libtorch, which the coverage CI job never is
# (it installs both on purpose; see test-coverage.yaml). Mocking
# torch_available()/requireNamespace() exercises both messages regardless of
# what is actually installed, so the branches are covered everywhere.
test_that("require_torch() reports the torch-package-missing message, mocked", {
  local_mocked_bindings(torch_available = function() FALSE)
  local_mocked_bindings(requireNamespace = function(package, ...) FALSE,
                         .package = "base")
  expect_error(require_torch(), "needs the 'torch' package")
  expect_error(require_torch(), "install\\.packages\\('torch'\\)")
  expect_error(require_torch(), "density_estimator = 'linear_gaussian'")
})

test_that("require_torch() reports the libtorch-missing message, mocked", {
  local_mocked_bindings(torch_available = function() FALSE)
  local_mocked_bindings(requireNamespace = function(package, ...) TRUE,
                         .package = "base")
  expect_error(require_torch(), "libtorch is not\\. Run torch::install_torch")
})

test_that("a neural estimator without torch fails with the install message", {
  # The user-visible route into require_torch(): check_torch_for_estimator()
  # (R/check.R) raises this before prepare_simulations() ever calls the
  # simulator (#250) -- asserted with a counting simulator below, not just
  # skipped past.
  skip_if(torch_available(), "torch works here, so the failing branch is unreachable")
  env <- new.env(parent = emptyenv())
  env$calls <- 0L
  sim <- function(theta) {
    env$calls <- env$calls + 1L
    theta + stats::rnorm(1, sd = 0.3)
  }
  expect_error(
    npe(prior_normal(mean = 0, sd = 1), sim,
        n_simulations = 100, density_estimator = "mdn", max_epochs = 2L),
    "'torch'"
  )
  expect_identical(env$calls, 0L)
})

test_that("npe()/nle()/nre() ask for torch before simulating, not after (#250)", {
  # Mocked, so this runs regardless of whether the test machine actually has
  # torch (the coverage CI job installs it on purpose; see test-coverage.yaml,
  # and test-utils.R's require_torch() tests above do the same). Before the
  # fix, all three burned the whole n_simulations budget in
  # prepare_simulations() before require_torch() was ever reached -- it lived
  # deep inside train_restarts(), which prepare_simulations() runs before.
  local_mocked_bindings(torch_available = function() FALSE)
  local_mocked_bindings(requireNamespace = function(package, ...) FALSE,
                        .package = "base")

  env <- new.env(parent = emptyenv())
  env$calls <- 0L
  sim <- function(theta) {
    env$calls <- env$calls + 1L
    theta + stats::rnorm(1, sd = 0.3)
  }
  prior <- prior_normal(mean = 0, sd = 1)

  # Defaults: npe()/nle() default to "maf", nre() to "resnet", all torch-only.
  expect_error(npe(prior, sim, n_simulations = 500),
               "needs the 'torch' package")
  expect_identical(env$calls, 0L)

  expect_error(nle(prior, sim, n_simulations = 500),
               "needs the 'torch' package")
  expect_identical(env$calls, 0L)

  expect_error(nre(prior, sim, n_simulations = 500),
               "needs the 'torch' package")
  expect_identical(env$calls, 0L)
  # nre()'s message points at its own torch-free alternative, not npe()'s.
  expect_error(nre(prior, sim, n_simulations = 500),
               'classifier = "logistic"')

  # The torch-free paths are untouched: they still simulate and fit.
  fit <- npe(prior, sim, n_simulations = 20,
            density_estimator = "linear_gaussian")
  expect_identical(fit$density_estimator, "linear_gaussian")
  expect_identical(env$calls, 20L)

  re <- nre(prior, sim, n_simulations = 20, classifier = "logistic")
  expect_identical(re$classifier, "logistic")
  expect_identical(env$calls, 40L)
})

test_that("npe()/nle()/nre() check device availability before simulating, not after (#250)", {
  # resolve_device()'s CUDA/MPS check needs torch loaded to even ask the
  # question, so this one is skipped rather than mocked when torch itself is
  # absent -- there would be nothing to mock cuda_is_available() onto.
  skip_if_no_torch()
  local_mocked_bindings(cuda_is_available = function() FALSE, .package = "torch")

  env <- new.env(parent = emptyenv())
  env$calls <- 0L
  sim <- function(theta) {
    env$calls <- env$calls + 1L
    theta + stats::rnorm(1, sd = 0.3)
  }
  prior <- prior_normal(mean = 0, sd = 1)

  expect_error(
    npe(prior, sim, n_simulations = 500, density_estimator = "mdn",
        device = "cuda", max_epochs = 2L),
    "cuda_is_available"
  )
  expect_identical(env$calls, 0L)

  # `device` is ignored by linear_gaussian, so an unavailable CUDA build is
  # never noticed -- and the simulator still runs.
  fit <- npe(prior, sim, n_simulations = 10,
            density_estimator = "linear_gaussian", device = "cuda")
  expect_identical(fit$density_estimator, "linear_gaussian")
  expect_identical(env$calls, 10L)
})

test_that("require_ggplot2() names the plotting package that is missing", {
  skip_if(requireNamespace("ggplot2", quietly = TRUE),
          "ggplot2 is installed, so the failing branch is unreachable")
  expect_error(require_ggplot2(), "needs the 'ggplot2' package")
  # ggally/ggdensity are only reached once ggplot2 itself is there, so with
  # ggplot2 missing the message is still about ggplot2.
  expect_error(require_ggplot2(ggally = TRUE, ggdensity = TRUE),
               "needs the 'ggplot2' package")
})

test_that("require_ggplot2() passes when the plotting stack is installed", {
  skip_if_no_ggplot2()
  expect_true(require_ggplot2())
  skip_if_no_ggally()
  expect_true(require_ggplot2(ggally = TRUE, ggdensity = TRUE))
})

# As with require_torch() above, the missing-GGally/ggdensity branches only
# run when those packages actually are missing, which the coverage CI job
# never is. Mock requireNamespace() per package so each branch is covered
# regardless of what is actually installed.
test_that("require_ggplot2() reports which optional package is missing, mocked", {
  local_mocked_bindings(requireNamespace = function(package, ...) !identical(package, "GGally"),
                         .package = "base")
  expect_error(require_ggplot2(ggally = TRUE), "needs the 'GGally' package")

  local_mocked_bindings(requireNamespace = function(package, ...) !identical(package, "ggdensity"),
                         .package = "base")
  expect_error(require_ggplot2(ggdensity = TRUE), "needs the 'ggdensity' package")
})

test_that("as_theta_matrix() coerces a data frame and keeps its column names", {
  df <- data.frame(alpha = c(1, 2), beta = c(3, 4))
  m <- as_theta_matrix(df, 2L)
  expect_true(is.matrix(m))
  expect_equal(unname(m), matrix(c(1, 2, 3, 4), nrow = 2L))
  expect_equal(colnames(m), c("alpha", "beta"))
})

test_that("verbose_cat() prints only when verbose is TRUE", {
  expect_silent(verbose_cat(FALSE, "hello\n"))
  expect_output(verbose_cat(TRUE, "hello\n"), "hello")
})
