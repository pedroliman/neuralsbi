#' Benchmark tasks
#'
#' Standard simulation-based-inference benchmark tasks, following the
#' definitions in the `sbibm` benchmark suite. Each task bundles a prior, a
#' simulator, and (where one exists) an analytic reference posterior, so the
#' same object drives unit tests, calibration studies, and the benchmark
#' harness in `inst/benchmarks/`.
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
#' @param dim Parameter/data dimension (default 10).
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

#' SIR epidemic task (sbibm `sir`)
#'
#' Susceptible-Infected-Recovered dynamics with contact rate \eqn{\beta} and
#' recovery rate \eqn{\gamma} under log-normal priors (as in `sbibm`):
#' \eqn{\beta \sim \mathrm{LogNormal}(\log 0.4, 0.5)},
#' \eqn{\gamma \sim \mathrm{LogNormal}(\log 0.125, 0.2)}. The ODE is solved
#' by simple Euler steps for a population of `N` over `days` days, and the
#' data are binomial subsamples (`n_obs_draws` trials) of the infected
#' fraction at `n_points` evenly spaced times. No closed-form posterior:
#' verify with SBC / expected coverage (Level 2).
#'
#' @param N,days,n_points,n_obs_draws Population size, horizon, number of
#'   observation times, and binomial trials per observation.
#' @rdname tasks
#' @export
task_sir <- function(N = 1e6, days = 160, n_points = 10L, n_obs_draws = 1000L) {
  prior <- prior_custom(
    sample_fn = function(n) cbind(stats::rlnorm(n, log(0.4), 0.5),
                                  stats::rlnorm(n, log(0.125), 0.2)),
    log_prob_fn = function(theta) {
      theta <- as_theta_matrix(theta, 2L)
      stats::dlnorm(theta[, 1], log(0.4), 0.5, log = TRUE) +
        stats::dlnorm(theta[, 2], log(0.125), 0.2, log = TRUE)
    },
    dim = 2L, lower = c(0, 0)
  )
  obs_times <- round(seq(1, days, length.out = n_points))
  simulator <- function(theta) {
    theta <- as_theta_matrix(theta, 2L)
    n <- nrow(theta)
    x <- matrix(0, nrow = n, ncol = n_points)
    for (i in seq_len(n)) {
      beta <- theta[i, 1]; gamma <- theta[i, 2]
      S <- N - 1; I <- 1; R <- 0
      Ipath <- numeric(days)
      for (t in seq_len(days)) {
        newinf <- beta * S * I / N
        newrec <- gamma * I
        S <- S - newinf
        I <- I + newinf - newrec
        R <- R + newrec
        Ipath[t] <- I
      }
      p <- pmin(pmax(Ipath[obs_times] / N, 0), 1)
      x[i, ] <- stats::rbinom(n_points, n_obs_draws, p) / n_obs_draws
    }
    x
  }
  new_task("sir", prior, simulator, 2L, n_points)
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
