#' Benchmark tasks
#'
#' Standard simulation-based-inference benchmark tasks, matching the
#' definitions in the `sbibm` suite that the Python `sbi` package benchmarks
#' against. Each task bundles a prior, a simulator, and (where one exists) an
#' analytic reference posterior, so the same object drives unit tests,
#' calibration studies, and the head-to-head comparison with Python `sbi`
#' (`inst/benchmarks/`).
#'
#' * [task_gaussian_linear()] -- conjugate Gaussian; analytic posterior.
#' * [task_two_moons()] -- crescent-shaped, bimodal posterior.
#' * [task_slcp()] -- "simple likelihood, complex posterior"; 5 parameters,
#'   8-dimensional data, strongly non-Gaussian posterior.
#'
#' @return A list of class `nsbi_task` with elements `name`, `prior`,
#'   `simulator`, `dim_theta`, `dim_x`, and optionally
#'   `reference_posterior(x_obs, n)` returning exact posterior draws.
#' @name tasks
NULL

#' @keywords internal
new_task <- function(name, prior, simulator, dim_theta, dim_x,
                     reference_posterior = NULL) {
  structure(
    list(name = name, prior = prior, simulator = simulator,
         dim_theta = dim_theta, dim_x = dim_x,
         reference_posterior = reference_posterior),
    class = "nsbi_task"
  )
}

#' @export
print.nsbi_task <- function(x, ...) {
  cat(sprintf("<nsbi_task> %s: %d parameters -> %d data dims%s\n",
              x$name, x$dim_theta, x$dim_x,
              if (is.null(x$reference_posterior)) "" else
                " (analytic reference available)"))
  invisible(x)
}

#' Gaussian linear task (sbibm `gaussian_linear`)
#'
#' Prior \eqn{\theta \sim N(0, 0.1 I)}, likelihood
#' \eqn{x \mid \theta \sim N(\theta, 0.1 I)}. The posterior is conjugate
#' Gaussian, so `reference_posterior()` is exact.
#'
#' @param dim Parameter/data dimension (sbibm uses 10).
#' @param prior_var,noise_var Prior and likelihood variances.
#' @rdname tasks
#' @export
task_gaussian_linear <- function(dim = 10L, prior_var = 0.1, noise_var = 0.1) {
  prior <- prior_normal(mean = rep(0, dim), sd = sqrt(prior_var))
  simulator <- function(theta) {
    theta <- as_theta_matrix(theta, dim)
    theta + matrix(stats::rnorm(length(theta), sd = sqrt(noise_var)),
                   nrow = nrow(theta))
  }
  reference <- function(x_obs, n = 10000L) {
    x_obs <- as.numeric(x_obs)
    post_var <- 1 / (1 / prior_var + 1 / noise_var)
    post_mean <- post_var * x_obs / noise_var
    matrix(stats::rnorm(n * dim, mean = rep(post_mean, each = n),
                        sd = sqrt(post_var)), nrow = n)
  }
  new_task("gaussian_linear", prior, simulator, dim, dim, reference)
}

#' Two-moons task (sbibm `two_moons`)
#'
#' Uniform prior on \eqn{[-1, 1]^2}; the simulator maps parameters onto a
#' noisy crescent whose posterior (for typical observations) has two
#' symmetric modes -- the standard bimodality stress test.
#'
#' @rdname tasks
#' @export
task_two_moons <- function() {
  prior <- prior_uniform(low = c(-1, -1), high = c(1, 1))
  simulator <- function(theta) {
    theta <- as_theta_matrix(theta, 2L)
    n <- nrow(theta)
    a <- stats::runif(n, -pi / 2, pi / 2)
    r <- stats::rnorm(n, mean = 0.1, sd = 0.01)
    px <- r * cos(a) + 0.25
    py <- r * sin(a)
    cbind(px - abs(theta[, 1] + theta[, 2]) / sqrt(2),
          py + (-theta[, 1] + theta[, 2]) / sqrt(2))
  }
  new_task("two_moons", prior, simulator, 2L, 2L)
}

#' SLCP task (sbibm `slcp`)
#'
#' "Simple likelihood, complex posterior": uniform prior on
#' \eqn{[-3, 3]^5}; the data are four i.i.d. draws from a 2-D Gaussian whose
#' mean, scales, and correlation are nonlinear functions of \eqn{\theta},
#' flattened to 8 dimensions. The posterior is multi-modal with heavy
#' nonlinear structure -- flows (MAF/NSF) are needed to fit it well.
#'
#' @rdname tasks
#' @export
task_slcp <- function() {
  prior <- prior_uniform(low = rep(-3, 5), high = rep(3, 5))
  simulator <- function(theta) {
    theta <- as_theta_matrix(theta, 5L)
    n <- nrow(theta)
    m1 <- theta[, 1]; m2 <- theta[, 2]
    s1 <- theta[, 3]^2; s2 <- theta[, 4]^2
    rho <- tanh(theta[, 5])
    x <- matrix(0, nrow = n, ncol = 8L)
    for (j in 0:3) {
      z1 <- stats::rnorm(n); z2 <- stats::rnorm(n)
      x[, 2 * j + 1L] <- m1 + s1 * z1
      x[, 2 * j + 2L] <- m2 + s2 * (rho * z1 + sqrt(1 - rho^2) * z2)
    }
    x
  }
  new_task("slcp", prior, simulator, 5L, 8L)
}
