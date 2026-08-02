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

test_that("a neural estimator without torch fails with the install message", {
  # The user-visible route into require_torch(): the controls are checked and
  # the simulations run first, then training asks for torch.
  skip_if(torch_available(), "torch works here, so the failing branch is unreachable")
  expect_error(
    npe(prior_normal(mean = 0, sd = 1),
        function(theta) theta + stats::rnorm(1, sd = 0.3),
        n_simulations = 100, density_estimator = "mdn", max_epochs = 2L),
    "'torch'"
  )
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
