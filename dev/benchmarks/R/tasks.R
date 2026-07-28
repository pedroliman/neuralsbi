# The ten sbibm benchmark tasks, reimplemented in R.
#
# Each task mirrors `sbibm/tasks/<name>/task.py` line for line: same prior, same
# simulator, same summary statistics, same data dimension. The simulators here
# are vectorized over parameter rows because the benchmark needs up to 100k
# draws; the one-row-at-a-time contract that `npe()` expects from a user
# simulator would make that unbearably slow, so the scripts pass `theta` and `x`
# to `npe()`/`nle()` directly.
#
# Two deliberate departures from sbibm, both documented in the README:
#
#   * SIR and Lotka-Volterra integrate their ODEs with a fixed-step vectorized
#     RK4 instead of Julia's adaptive Tsit5. The step sizes below put the
#     integration error orders of magnitude below the observation noise.
#   * Everything else is exact, including the frozen constants (GLM stimulus,
#     SLCP distractor noise) that are read from sbibm's own `.pt` files.

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

#' One vectorized RK4 step over a list-of-vectors state.
rk4_step <- function(state, deriv, h) {
  k1 <- deriv(state)
  s2 <- Map(function(u, k) u + 0.5 * h * k, state, k1)
  k2 <- deriv(s2)
  s3 <- Map(function(u, k) u + 0.5 * h * k, state, k2)
  k3 <- deriv(s3)
  s4 <- Map(function(u, k) u + h * k, state, k3)
  k4 <- deriv(s4)
  Map(function(u, a, b, c, d) u + (h / 6) * (a + 2 * b + 2 * c + d),
      state, k1, k2, k3, k4)
}

#' Integrate `state` forward, recording it at each of `n_out` equally spaced
#' output times separated by `interval`.
integrate_grid <- function(state, deriv, interval, n_out, steps_per_interval) {
  h <- interval / steps_per_interval
  out <- vector("list", n_out)
  out[[1]] <- state
  for (j in seq_len(n_out - 1)) {
    for (s in seq_len(steps_per_interval)) state <- rk4_step(state, deriv, h)
    out[[j + 1]] <- state
  }
  out
}

task_sir_ <- function(N = 1e6, I0 = 1, R0 = 0, steps_per_interval = 170L) {
  prior <- prior_lognormal(c(log(0.4), log(0.125)), c(0.5, 0.2))
  simulate <- function(theta) {
    n <- nrow(theta)
    beta <- theta[, 1]; gamma <- theta[, 2]
    deriv <- function(u) {
      inf <- beta * u$S * u$I / N
      rec <- gamma * u$I
      list(S = -inf, I = inf - rec, R = rec)
    }
    # sbibm saves daily and then keeps every 17th day: 10 points from day 0.
    grid <- integrate_grid(list(S = rep(N - I0 - R0, n), I = rep(I0, n),
                                R = rep(R0, n)),
                           deriv, interval = 17, n_out = 10L,
                           steps_per_interval = steps_per_interval)
    I <- vapply(grid, function(u) u$I, numeric(n))
    if (n == 1) I <- matrix(I, nrow = 1)
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

task_lv <- function(steps_per_interval = 200L) {
  prior <- prior_lognormal(c(-0.125, -3, -0.125, -3), rep(0.5, 4))
  simulate <- function(theta) {
    n <- nrow(theta)
    alpha <- theta[, 1]; beta <- theta[, 2]
    gamma <- theta[, 3]; delta <- theta[, 4]
    deriv <- function(u) {
      list(x = alpha * u$x - beta * u$x * u$y,
           y = -gamma * u$y + delta * u$x * u$y)
    }
    # saveat = 0.1 over 20 days, then every 21st point: 10 points, 2.1 apart.
    grid <- integrate_grid(list(x = rep(30, n), y = rep(1, n)),
                           deriv, interval = 2.1, n_out = 10L,
                           steps_per_interval = steps_per_interval)
    xs <- vapply(grid, function(u) u$x, numeric(n))
    ys <- vapply(grid, function(u) u$y, numeric(n))
    if (n == 1) { xs <- matrix(xs, nrow = 1); ys <- matrix(ys, nrow = 1) }
    us <- cbind(xs, ys)                     # sbibm flattens species-major
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
