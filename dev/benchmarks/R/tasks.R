# The ten sbibm benchmark tasks, reimplemented in R.
#
# Each task mirrors `sbibm/tasks/<name>/task.py` line for line: same prior, same
# simulator, same summary statistics, same data dimension. The simulators here
# are vectorized over parameter rows because the benchmark needs up to 100k
# draws; the one-row-at-a-time contract that `npe()` expects from a user
# simulator would make that unbearably slow, so the scripts pass `theta` and `x`
# to `npe()`/`nle()` directly.
#
# SIR and Lotka-Volterra are the exception to the vectorization: sbibm solves
# their ODEs one parameter set at a time with Julia's DifferentialEquations, and
# we do the same with deSolve's `lsoda`, which is the closest thing R has to
# that solver (adaptive step, automatic stiff/non-stiff switching). The noise
# model on top of the trajectories is still vectorized.
#
# Everything else is exact, including the frozen constants (GLM stimulus, SLCP
# distractor noise) that are read from sbibm's own `.pt` files.

sbibm_task_names <- function() {
  c("gaussian_linear", "gaussian_linear_uniform", "gaussian_mixture",
    "two_moons", "slcp", "slcp_distractors", "bernoulli_glm",
    "bernoulli_glm_raw", "sir", "lotka_volterra")
}

#' Build one of the sbibm tasks.
#'
#' @return A list with `name`, `dim_theta`, `dim_x`, an `nsbi_prior`, a
#'   vectorized `simulate(theta)`, and the sbibm file locations for the
#'   observations and reference posteriors.
sbibm_task <- function(name) {
  name <- match.arg(name, sbibm_task_names())
  switch(name,
    gaussian_linear         = task_gl(),
    gaussian_linear_uniform = task_glu(),
    gaussian_mixture        = task_gm(),
    two_moons               = task_tm(),
    slcp                    = task_slcp_(FALSE),
    slcp_distractors        = task_slcp_(TRUE),
    bernoulli_glm           = task_glm(FALSE),
    bernoulli_glm_raw       = task_glm(TRUE),
    sir                     = task_sir_(),
    lotka_volterra          = task_lv()
  )
}

new_bench_task <- function(name, dim_theta, dim_x, prior, simulate, files_dir,
                           observation_file = "observation.csv",
                           budgets = c(1e3, 1e4, 1e5)) {
  structure(list(name = name, dim_theta = dim_theta, dim_x = dim_x,
                 prior = prior, simulate = simulate, files_dir = files_dir,
                 observation_file = observation_file, budgets = budgets),
            class = "sbibm_bench_task")
}

print.sbibm_bench_task <- function(x, ...) {
  cat(sprintf("<sbibm task> %s: %d parameters -> %d data dims\n",
              x$name, x$dim_theta, x$dim_x))
  invisible(x)
}

# --- priors sbibm needs that neuralsbi does not ship -------------------------

#' Multivariate normal prior given a covariance matrix.
prior_mvnorm <- function(mean, cov, lower = NULL, upper = NULL) {
  d <- length(mean)
  L <- chol(cov)                       # upper triangular, cov = t(L) %*% L
  prec <- chol2inv(L)
  logdet <- 2 * sum(log(diag(L)))
  prior_custom(
    sample_fn = function(n) {
      z <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
      sweep(z %*% L, 2, mean, `+`)
    },
    log_prob_fn = function(theta) {
      theta <- as.matrix(theta)
      if (ncol(theta) != d) theta <- matrix(theta, ncol = d, byrow = TRUE)
      dev <- sweep(theta, 2, mean, `-`)
      -0.5 * (rowSums((dev %*% prec) * dev) + logdet + d * log(2 * pi))
    },
    dim = d, lower = lower, upper = upper
  )
}

#' Independent log-normal prior.
prior_lognormal <- function(meanlog, sdlog) {
  d <- length(meanlog)
  prior_custom(
    sample_fn = function(n) {
      matrix(stats::rlnorm(n * d, rep(meanlog, each = n), rep(sdlog, each = n)),
             nrow = n, ncol = d)
    },
    log_prob_fn = function(theta) {
      theta <- as.matrix(theta)
      if (ncol(theta) != d) theta <- matrix(theta, ncol = d, byrow = TRUE)
      lp <- matrix(stats::dlnorm(theta, rep(meanlog, each = nrow(theta)),
                                 rep(sdlog, each = nrow(theta)), log = TRUE),
                   nrow = nrow(theta))
      out <- rowSums(lp)
      out[!is.finite(out)] <- -Inf
      out
    },
    dim = d, lower = rep(0, d), upper = rep(Inf, d)
  )
}

# --- gaussian_linear ---------------------------------------------------------

# NOTE the parameterisation: sbibm passes `prior_scale` and `simulator_scale`
# through `torch.inverse(scale * eye)`, so 0.1 is a *variance*, not an sd.
task_gl <- function(dim = 10L, prior_var = 0.1, noise_var = 0.1) {
  prior <- prior_normal(mean = rep(0, dim), sd = sqrt(prior_var))
  simulate <- function(theta) {
    theta + matrix(stats::rnorm(length(theta), sd = sqrt(noise_var)),
                   nrow = nrow(theta), ncol = dim)
  }
  new_bench_task("gaussian_linear", dim, dim, prior, simulate, "gaussian_linear")
}

task_glu <- function(dim = 10L, prior_bound = 1, noise_var = 0.1) {
  prior <- prior_uniform(rep(-prior_bound, dim), rep(prior_bound, dim))
  simulate <- function(theta) {
    theta + matrix(stats::rnorm(length(theta), sd = sqrt(noise_var)),
                   nrow = nrow(theta), ncol = dim)
  }
  new_bench_task("gaussian_linear_uniform", dim, dim, prior, simulate,
                 "gaussian_linear_uniform")
}

# --- gaussian_mixture --------------------------------------------------------

task_gm <- function(prior_bound = 10) {
  prior <- prior_uniform(rep(-prior_bound, 2), rep(prior_bound, 2))
  scales <- c(1, 0.1)
  simulate <- function(theta) {
    n <- nrow(theta)
    # One mixture component per parameter row; both coordinates share it.
    idx <- sample.int(2L, n, replace = TRUE)
    sd <- scales[idx]
    theta + matrix(stats::rnorm(2 * n), nrow = n, ncol = 2) * sd
  }
  new_bench_task("gaussian_mixture", 2L, 2L, prior, simulate, "gaussian_mixture")
}

# --- two_moons ---------------------------------------------------------------

task_tm <- function() {
  prior <- prior_uniform(c(-1, -1), c(1, 1))
  ang <- -pi / 4
  cc <- cos(ang); ss <- sin(ang)
  simulate <- function(theta) {
    n <- nrow(theta)
    a <- stats::runif(n, -pi / 2, pi / 2)
    r <- stats::rnorm(n, 0.1, 0.01)
    p1 <- cos(a) * r + 0.25
    p2 <- sin(a) * r
    z0 <- cc * theta[, 1] - ss * theta[, 2]
    z1 <- ss * theta[, 1] + cc * theta[, 2]
    cbind(p1 - abs(z0), p2 + z1)
  }
  new_bench_task("two_moons", 2L, 2L, prior, simulate, "two_moons")
}

# --- slcp / slcp_distractors -------------------------------------------------

task_slcp_ <- function(distractors) {
  prior <- prior_uniform(rep(-3, 5), rep(3, 5))
  eps <- 1e-6
  sim_core <- function(theta) {
    n <- nrow(theta)
    s1 <- theta[, 3]^2
    s2 <- theta[, 4]^2
    rho <- tanh(theta[, 5])
    a <- s1^2 + eps                      # var of coordinate 1
    d <- s2^2 + eps                      # var of coordinate 2
    b <- rho * s1 * s2                   # covariance
    # Cholesky of the 2x2 covariance, done in closed form.
    l11 <- sqrt(a)
    l21 <- b / l11
    l22 <- sqrt(pmax(d - l21^2, 0))
    out <- matrix(0, nrow = n, ncol = 8)
    for (k in 0:3) {                     # four i.i.d. draws per parameter set
      z1 <- stats::rnorm(n)
      z2 <- stats::rnorm(n)
      out[, 2 * k + 1] <- theta[, 1] + l11 * z1
      out[, 2 * k + 2] <- theta[, 2] + l21 * z1 + l22 * z2
    }
    out
  }
  if (!distractors) {
    return(new_bench_task("slcp", 5L, 8L, prior, sim_core, "slcp"))
  }
  simulate <- function(theta) {
    x <- sim_core(theta)
    noise <- slcp_noise_sample(nrow(theta))
    cbind(x, noise)[, slcp_permutation(), drop = FALSE]
  }
  new_bench_task("slcp_distractors", 5L, 100L, prior, simulate, "slcp",
                 observation_file = "observation_distractors.csv")
}

.slcp_cache <- new.env(parent = emptyenv())

#' Parameters of the 20-component multivariate-t mixture sbibm adds as
#' distractor dimensions, read from its `gmm.torch`.
slcp_noise_params <- function() {
  if (!is.null(.slcp_cache$gmm)) return(.slcp_cache$gmm)
  path <- file.path(sbibm_root(), "sbibm", "tasks", "slcp", "files", "gmm.torch")
  loc <- matrix(read_torch_storage(path, 20 * 92, "double", which = "first"),
                nrow = 20, ncol = 92, byrow = TRUE)
  tril_flat <- read_torch_storage(path, 20 * 92 * 92, "double", which = "last")
  scale_tril <- lapply(seq_len(20), function(k) {
    off <- (k - 1) * 92 * 92
    matrix(tril_flat[(off + 1):(off + 92 * 92)], nrow = 92, ncol = 92,
           byrow = TRUE)
  })
  .slcp_cache$gmm <- list(loc = loc, scale_tril = scale_tril, df = 2)
  .slcp_cache$gmm
}

slcp_permutation <- function() {
  if (!is.null(.slcp_cache$perm)) return(.slcp_cache$perm)
  path <- file.path(sbibm_root(), "sbibm", "tasks", "slcp", "files",
                    "permutation_idx.torch")
  perm <- read_torch_storage(path, 100, "int64", which = "last") + 1
  stopifnot(setequal(perm, 1:100))
  .slcp_cache$perm <- perm
  perm
}

#' Draw n rows from the distractor mixture.
#'
#' pyro's MultivariateStudentT samples as loc + L %*% (z * sqrt(df / chi2_df),
#' with a uniform categorical over the 20 components.
slcp_noise_sample <- function(n) {
  g <- slcp_noise_params()
  idx <- sample.int(20L, n, replace = TRUE)
  out <- matrix(0, nrow = n, ncol = 92)
  scale <- sqrt(g$df / stats::rchisq(n, df = g$df))
  for (k in unique(idx)) {
    rows <- which(idx == k)
    z <- matrix(stats::rnorm(length(rows) * 92), nrow = length(rows), ncol = 92)
    y <- z * scale[rows]
    out[rows, ] <- sweep(y %*% t(g$scale_tril[[k]]), 2, g$loc[k, ], `+`)
  }
  out
}

# --- bernoulli_glm / bernoulli_glm_raw ---------------------------------------

.glm_cache <- new.env(parent = emptyenv())

#' sbibm's frozen 100 x 10 design matrix: a column of ones followed by nine
#' lagged copies of the Gaussian white-noise stimulus.
glm_design_matrix <- function() {
  if (!is.null(.glm_cache$X)) return(.glm_cache$X)
  path <- file.path(sbibm_root(), "sbibm", "tasks", "bernoulli_glm", "files",
                    "design_matrix.pt")
  X <- matrix(read_torch_storage(path, 1000, "float", which = "last"),
              nrow = 100, ncol = 10, byrow = TRUE)
  stopifnot(all(X[, 1] == 1))
  .glm_cache$X <- X
  X
}

glm_prior <- function() {
  M <- 9L
  D <- diag(M)
  D[cbind(2:M, 1:(M - 1))] <- -1
  Fm <- D %*% D + diag(sqrt((0:(M - 1)) / M))
  Binv <- matrix(0, M + 1, M + 1)
  Binv[1, 1] <- 0.5                        # offset
  Binv[-1, -1] <- t(Fm) %*% Fm             # smoothness penalty on the filter
  prior_mvnorm(rep(0, M + 1), solve(Binv))
}

task_glm <- function(raw) {
  prior <- glm_prior()
  simulate <- function(theta) {
    X <- glm_design_matrix()
    n <- nrow(theta)
    psi <- theta %*% t(X)                  # n x 100
    z <- 1 / (1 + exp(-psi))
    y <- matrix(as.numeric(stats::runif(n * 100) < z), nrow = n, ncol = 100)
    if (raw) return(y)
    # The sufficient statistics are the spike count and the nine-lag
    # spike-triggered average, which together are exactly `y %*% X`: column 1 of
    # the design matrix is the ones vector, columns 2-10 the lagged stimulus.
    y %*% X
  }
  if (raw) {
    new_bench_task("bernoulli_glm_raw", 10L, 100L, prior, simulate,
                   "bernoulli_glm", observation_file = "observation_raw.csv")
  } else {
    new_bench_task("bernoulli_glm", 10L, 10L, prior, simulate, "bernoulli_glm")
  }
}

# --- ODE tasks ---------------------------------------------------------------

require_desolve <- function() {
  if (!requireNamespace("deSolve", quietly = TRUE)) {
    stop("The sir and lotka_volterra tasks need deSolve: ",
         'install.packages("deSolve")', call. = FALSE)
  }
}

#' Solve one ODE per row of `theta` and return the trajectories stacked.
#'
#' sbibm solves these one parameter set at a time and replaces a trajectory it
#' could not integrate with NaN, which then propagates into the observation.
#' Same here: a solve that errors, stops short of the last requested time, or
#' returns anything non-finite becomes a row of NaN.
#'
#' @param func A deSolve derivative function `function(t, u, parms)`.
#' @param y0 Initial state.
#' @param times Output grid, matching sbibm's `saveat`.
#' @return An `n x (length(times) * length(y0))` matrix, state-major within each
#'   row (all times of state 1, then all times of state 2), which is how sbibm
#'   flattens `num_parameters x num_states x num_times`.
solve_ode_rows <- function(theta, func, y0, times) {
  require_desolve()
  n_state <- length(y0)
  n_time <- length(times)
  blank <- rep(NA_real_, n_state * n_time)

  one <- function(i) {
    out <- tryCatch(
      suppressWarnings(
        deSolve::ode(y = y0, times = times, func = func, parms = theta[i, ])
      ),
      error = function(e) NULL
    )
    if (is.null(out) || nrow(out) != n_time) return(blank)
    u <- out[, -1, drop = FALSE]           # drop the time column
    if (!all(is.finite(u))) return(blank)
    as.numeric(u)                          # column-major: state-major, as wanted
  }

  rows <- bench_lapply(seq_len(nrow(theta)), one)
  matrix(unlist(rows, use.names = FALSE), nrow = nrow(theta), byrow = TRUE)
}

#' `lapply`, or `mclapply` when NEURALSBI_BENCH_CORES asks for it.
#'
#' Only the deterministic ODE solves go through this. The stochastic part of
#' each simulator stays in the parent process, so the number of cores does not
#' change the draws for a given seed.
bench_lapply <- function(x, f) {
  cores <- suppressWarnings(as.integer(Sys.getenv("NEURALSBI_BENCH_CORES", "1")))
  if (is.na(cores) || cores <= 1 || .Platform$OS.type == "windows") {
    return(lapply(x, f))
  }
  parallel::mclapply(x, f, mc.cores = cores)
}

task_sir_ <- function(N = 1e6, I0 = 1, R0 = 0) {
  prior <- prior_lognormal(c(log(0.4), log(0.125)), c(0.5, 0.2))
  # saveat = 1 day over 160 days, then every 17th point: 10 points from day 0.
  times <- seq(0, 160, by = 1)
  keep <- seq(1, length(times), by = 17)[1:10]
  deriv <- function(t, u, p) {
    inf <- p[1] * u[1] * u[2] / N
    rec <- p[2] * u[2]
    list(c(-inf, inf - rec, rec))
  }
  simulate <- function(theta) {
    n <- nrow(theta)
    us <- solve_ode_rows(theta, deriv, c(N - I0 - R0, I0, R0), times)
    # Columns are state-major, so the infected compartment is the second block.
    I <- us[, length(times) + keep, drop = FALSE]
    bad <- !is.finite(rowSums(I))
    p <- pmin(pmax(I / N, 0), 1)
    p[bad, ] <- 0
    x <- matrix(stats::rbinom(length(p), size = 1000, prob = p),
                nrow = n, ncol = 10)
    x[bad, ] <- NA_real_
    x
  }
  new_bench_task("sir", 2L, 10L, prior, simulate, "sir")
}

task_lv <- function() {
  prior <- prior_lognormal(c(-0.125, -3, -0.125, -3), rep(0.5, 4))
  # saveat = 0.1 over 20 days, then every 21st point: 10 points, 2.1 apart.
  times <- seq(0, 20, by = 0.1)
  keep <- seq(1, length(times), by = 21)[1:10]
  deriv <- function(t, u, p) {
    list(c(p[1] * u[1] - p[2] * u[1] * u[2],
           -p[3] * u[2] + p[4] * u[1] * u[2]))
  }
  simulate <- function(theta) {
    n <- nrow(theta)
    all_states <- solve_ode_rows(theta, deriv, c(30, 1), times)
    us <- all_states[, c(keep, length(times) + keep), drop = FALSE]
    bad <- !is.finite(rowSums(us))
    clamped <- pmin(pmax(us, 1e-10), 1e4)
    clamped[bad, ] <- 1
    x <- exp(log(clamped) + matrix(stats::rnorm(length(clamped), sd = 0.1),
                                   nrow = n))
    x[bad, ] <- NA_real_
    x
  }
  new_bench_task("lotka_volterra", 4L, 20L, prior, simulate, "lotka_volterra")
}

# --- simulation driver -------------------------------------------------------

#' Run a task's simulator over `theta`, in chunks so peak memory stays bounded.
simulate_chunked <- function(task, theta, chunk = 20000L, verbose = FALSE) {
  n <- nrow(theta)
  if (n <= chunk) return(task$simulate(theta))
  starts <- seq(1, n, by = chunk)
  parts <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    rows <- starts[i]:min(starts[i] + chunk - 1, n)
    parts[[i]] <- task$simulate(theta[rows, , drop = FALSE])
    if (verbose) say(sprintf("  simulated %d / %d", max(rows), n))
  }
  do.call(rbind, parts)
}
