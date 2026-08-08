# device = "cpu"/"cuda"/"mps"/"gpu"/"auto" on npe()/nle() (#82). CPU is the
# default and needs no torch; anything else needs torch loaded to resolve, so
# most of these are torch-free checks of the syntactic validation and of the
# "linear_gaussian ignores device" no-op, with a smaller torch-gated group for
# the parts that actually need a device to resolve against.

toy_prior <- function() prior_uniform(c(mu = -2), c(mu = 2))
toy_simulator <- function(mu) c(y = mu + stats::rnorm(1, sd = 0.1))

counting_simulator <- function(env) {
  function(mu) {
    env$calls <- env$calls + 1L
    c(y = mu + stats::rnorm(1, sd = 0.1))
  }
}

test_that("check_device_arg() accepts the recognized keywords", {
  for (d in c("cpu", "cuda", "mps", "gpu", "auto")) {
    expect_identical(check_device_arg(d), d)
  }
})

test_that("check_device_arg() rejects anything else, naming the keywords", {
  expect_error(check_device_arg("gpu2"),
               "`device` must be one of \"cpu\", \"cuda\", \"mps\", \"gpu\", \"auto\"")
  expect_error(check_device_arg("CUDA"), "must be one of")
  expect_error(check_device_arg(c("cpu", "cuda")), "must be one of")
  expect_error(check_device_arg(1), "must be one of")
  expect_error(check_device_arg(NA_character_), "must be one of")
  expect_error(check_device_arg(NULL), "must be one of")
})

test_that("npe()/nle() reject a bad device before simulating", {
  env <- new.env(parent = emptyenv())
  env$calls <- 0L
  expect_error(npe(toy_prior(), counting_simulator(env), n_simulations = 100,
                   density_estimator = "linear_gaussian", device = "tpu"),
               "`device` must be one of")
  expect_identical(env$calls, 0L)

  env$calls <- 0L
  expect_error(nle(toy_prior(), counting_simulator(env), n_simulations = 100,
                   density_estimator = "linear_gaussian", device = "tpu"),
               "`device` must be one of")
  expect_identical(env$calls, 0L)
})

test_that("device defaults to \"cpu\" and needs no torch for linear_gaussian", {
  fit <- npe(toy_prior(), toy_simulator, n_simulations = 200,
             density_estimator = "linear_gaussian", seed = 1)
  expect_identical(fit$device, "cpu")
})

test_that("device is a no-op for linear_gaussian, not an error, even when unavailable", {
  # fit_linear_gaussian() has no `device` formal, so fit_density_estimator()'s
  # dots-intersection silently drops it (the same mechanism n_components/hidden
  # rely on) and check_device_arg() never touches torch. This has to hold even
  # on a machine with no CUDA/MPS at all, which this sandbox is one example of.
  fit <- npe(toy_prior(), toy_simulator, n_simulations = 200,
             density_estimator = "linear_gaussian", device = "cuda", seed = 1)
  expect_s3_class(fit, "nsbi_npe")
  expect_identical(fit$density_estimator, "linear_gaussian")

  fit_mps <- nle(toy_prior(), toy_simulator, n_simulations = 200,
                 density_estimator = "linear_gaussian", device = "mps", seed = 1)
  expect_s3_class(fit_mps, "nsbi_nle")
})

test_that("device is a no-op for a custom density_estimator function too", {
  custom <- function(theta, x) fit_linear_gaussian(theta, x)
  fit <- npe(toy_prior(), toy_simulator, n_simulations = 200,
             density_estimator = custom, device = "cuda", seed = 1)
  expect_s3_class(fit, "nsbi_npe")
  expect_identical(fit$density_estimator, "custom")
})

test_that("resolve_device() passes cpu straight through without touching torch", {
  # cpu never calls cuda_is_available()/backends_mps_is_available(), so this
  # branch is exercised whether or not torch is installed here.
  expect_identical(resolve_device("cpu"), "cpu")
})

test_that("device = \"cuda\"/\"mps\" error clearly when unavailable", {
  skip_if_no_torch()
  skip_if(torch::cuda_is_available(), "cuda is actually available here")
  expect_error(resolve_device("cuda"), "device = \"cuda\"")
  expect_error(resolve_device("cuda"), "cuda_is_available\\(\\) is FALSE")

  skip_if(torch::backends_mps_is_available(), "mps is actually available here")
  expect_error(resolve_device("mps"), "device = \"mps\"")
  expect_error(resolve_device("mps"), "backends_mps_is_available\\(\\) is FALSE")
})

test_that("device = \"gpu\"/\"auto\" resolve CUDA -> MPS -> CPU without erroring", {
  skip_if_no_torch()
  expected <- if (isTRUE(torch::cuda_is_available())) "cuda" else
    if (isTRUE(torch::backends_mps_is_available())) "mps" else "cpu"
  expect_identical(resolve_device("gpu"), expected)
  expect_identical(resolve_device("auto"), expected)
})

test_that("npe(..., density_estimator = 'mdn', device = 'cuda') errors clearly here", {
  skip_if_no_torch()
  skip_if(torch::cuda_is_available(), "cuda is actually available here")
  expect_error(
    npe(toy_prior(), toy_simulator, n_simulations = 100,
        density_estimator = "mdn", device = "cuda", max_epochs = 2L),
    "device = \"cuda\"")
})

test_that("npe(..., device = 'cpu') trains and records the device it used", {
  skip_if_no_torch()
  set.seed(5)
  fit <- npe(toy_prior(), toy_simulator, n_simulations = 300,
             density_estimator = "mdn", n_components = 1L, hidden = c(8L),
             max_epochs = 5L, device = "cpu", seed = 5)
  expect_identical(fit$device, "cpu")
  expect_identical(fit$de$device, "cpu")
  # torch_device objects wrap an external pointer, so two separately
  # constructed "cpu" devices are never identical() even though they name the
  # same device; torch_device has its own `==` for exactly this comparison.
  expect_true(net_device(fit$de$net) == torch::torch_device("cpu"))

  # de_log_prob()/de_sample() still work end to end (the ordinary path every
  # other neural test exercises, now routed through net_device()).
  lp <- de_log_prob(fit$de, 0, 0.1)
  expect_true(is.finite(lp))
  draws <- de_sample(fit$de, 0.1, 20)
  expect_equal(dim(draws), c(20L, 1L))
})

test_that("train_conditional_de() checks its controls before device resolution", {
  # check_train_controls() needs no torch; it has to fire before
  # resolve_device() does, so a bad batch_size is reported even when the
  # device asked for would itself have errored.
  theta <- matrix(stats::rnorm(50), ncol = 1)
  x <- matrix(stats::rnorm(50), ncol = 1)
  expect_error(
    train_conditional_de(build_net = function() stop("not reached"),
                         log_prob_fn = function(...) stop("not reached"),
                         theta = theta, x = x, batch_size = 0,
                         device = "not-a-real-device"),
    "`batch_size` must be")
})

test_that("load_npe() lands a reloaded fit's device field back on cpu", {
  skip_if_no_torch()
  set.seed(6)
  fit <- npe(toy_prior(), toy_simulator, n_simulations = 200,
             density_estimator = "mdn", n_components = 1L, hidden = c(8L),
             max_epochs = 3L, device = "cpu", seed = 6)
  # Simulate what a device = "mps"/"cuda" fit's recorded device would say,
  # since this sandbox cannot actually produce one to round-trip.
  fit$de$device <- "mps"
  fit$device <- "mps"

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  save_npe(fit, path)
  fit2 <- load_npe(path)

  expect_identical(fit2$de$device, "cpu")
  expect_identical(fit2$device, "cpu")
  expect_true(net_device(fit2$de$net) == torch::torch_device("cpu"))
})
