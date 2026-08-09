# Argument guards on the plotting functions. Drawing needs ggplot2, and
# pairplot() also needs GGally and ggdensity; all three are Suggests, so these
# use the skips from helper-suggests.R.

sbc_fixture <- function(seed = 1) {
  set.seed(seed)
  prior <- prior_normal(mean = c(mu = 0, nu = 0), sd = 1)
  simulator <- function(mu, nu) c(a = mu, b = nu) + stats::rnorm(2, sd = 0.4)
  fit <- npe(prior, simulator, n_simulations = 600,
             density_estimator = "linear_gaussian")
  # 100 trials over the 20 rank bins keeps chisq.test() off its
  # small-expected-count warning.
  sbc(fit, simulator, n_sbc = 100L, n_posterior_samples = 50L, seed = seed)
}

test_that("plot_sbc() refuses a parameter it cannot subscript", {
  skip_if_no_ggplot2()
  res <- sbc_fixture()

  expect_error(plot_sbc(res, param = 99),
               regexp = "`param` must be one parameter index between 1 and 2")
  # A name is an accepted value, so the names belong in the message.
  expect_error(plot_sbc(res, param = 99), regexp = "or one of mu, nu")
  expect_error(plot_sbc(res, param = 0), regexp = "between 1 and 2")
  expect_error(plot_sbc(res, param = 1.5), regexp = "between 1 and 2")
  expect_error(plot_sbc(res, param = c(1, 2)),
               regexp = "not a length-2 numeric vector")
  expect_error(plot_sbc(res, param = NA), regexp = "not NA")
})

test_that("plot_sbc() resolves a parameter name against the rank columns", {
  skip_if_no_ggplot2()
  res <- sbc_fixture()

  expect_error(plot_sbc(res, param = "sigma"),
               regexp = "not one of the parameter names: mu, nu")

  unnamed <- res
  colnames(unnamed$ranks) <- NULL
  expect_error(plot_sbc(unnamed, param = "mu"),
               regexp = "the parameters are unnamed")

  path <- tempfile(fileext = ".png")
  grDevices::png(path)
  by_name <- plot_sbc(res, param = "nu")
  by_index <- plot_sbc(res, param = 2L)
  grDevices::dev.off()

  expect_identical(by_name, by_index)
  expect_length(by_name, 100L)
  expect_true(file.exists(path))
  unlink(path)
})

test_that("pairplot() checks limits against the number of parameters", {
  skip_if_no_ggally()
  set.seed(2)
  draws <- matrix(stats::rnorm(1500), ncol = 3)

  expect_error(pairplot(draws, limits = list(c(-3, 3), c(-3, 3))),
               regexp = "`limits` has 2 elements but `samples` has 3 parameters")
  expect_error(pairplot(draws, limits = matrix(c(-3, 3), nrow = 1)),
               regexp = "`limits` has 1 row but `samples` has 3 parameters")
  expect_error(pairplot(draws, limits = c(-3, 3)),
               regexp = "`limits` must be a list of 3 elements")
})

test_that("pairplot() accepts one limit pair per parameter", {
  skip_if_no_ggally()
  set.seed(3)
  draws <- matrix(stats::rnorm(1500), ncol = 3)
  path <- tempfile(fileext = ".png")

  grDevices::png(path)
  from_list <- pairplot(draws, limits = list(c(-3, 3), c(-3, 3), c(-3, 3)))
  from_matrix <- pairplot(draws, limits = matrix(rep(c(-3, 3), each = 3),
                                                 nrow = 3))
  grDevices::dev.off()

  expect_s3_class(from_list, "ggmatrix")
  expect_s3_class(from_matrix, "ggmatrix")
  expect_true(file.exists(path))
  unlink(path)
})

panel_range <- function(ggmatrix, row, col, axis = c("x", "y")) {
  axis <- match.arg(axis)
  built <- ggplot2::ggplot_build(ggmatrix[row, col])
  built$layout$panel_params[[1]][[paste0(axis, ".range")]]
}

test_that("pairplot() aligns each parameter's axis range across its panels by default", {
  skip_if_no_ggally()
  set.seed(4)
  draws <- matrix(stats::rnorm(3000), ncol = 3)
  path <- tempfile(fileext = ".png")

  grDevices::png(path)
  p <- pairplot(draws)
  grDevices::dev.off()
  unlink(path)

  # Column 1 (theta[1]): the two lower-triangle panels below the diagonal and
  # the diagonal panel itself all plot theta[1] on the x axis, and used to
  # each pick their own range because geom_hdr()/geom_density() estimate
  # their density grid per panel.
  col1_x <- list(
    panel_range(p, 2, 1, "x"),
    panel_range(p, 3, 1, "x"),
    panel_range(p, 1, 1, "x")
  )
  expect_equal(col1_x[[2]], col1_x[[1]])
  expect_equal(col1_x[[3]], col1_x[[1]])

  # Row 3 (theta[3]): the two lower-triangle panels in that row and the
  # diagonal panel all plot theta[3] on the y (off-diagonal) or x (diagonal)
  # axis and must agree with each other.
  row3 <- list(
    panel_range(p, 3, 1, "y"),
    panel_range(p, 3, 2, "y"),
    panel_range(p, 3, 3, "x")
  )
  expect_equal(row3[[2]], row3[[1]])
  expect_equal(row3[[3]], row3[[1]])

  # Column 2 (theta[2]): panel[2,1]'s y axis and panel[3,2]'s x axis both
  # plot theta[2] and must agree.
  expect_equal(panel_range(p, 3, 2, "x"), panel_range(p, 2, 1, "y"))
})

test_that("pairplot() draws truth markers as a subdued grey, not a saturated red", {
  skip_if_no_ggally()
  set.seed(5)
  draws <- matrix(stats::rnorm(600), ncol = 2)
  path <- tempfile(fileext = ".png")

  grDevices::png(path)
  p <- pairplot(draws, truth = c(0, 0))
  grDevices::dev.off()
  unlink(path)

  built <- ggplot2::ggplot_build(p[2, 1])
  vline_colour <- built$data[[2]]$colour
  point_colour <- built$data[[3]]$colour
  expect_true(all(vline_colour == "grey30"))
  expect_true(all(point_colour == "grey30"))
  expect_false(any(c(vline_colour, point_colour) == "firebrick"))
})

test_that("pad_range() pads a non-constant range by 5% and a constant one by a fixed amount", {
  expect_equal(pad_range(c(0, 10)), c(-0.5, 10.5))
  expect_equal(pad_range(rep(3, 5)), c(2.5, 3.5))
})
