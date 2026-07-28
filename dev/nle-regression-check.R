## Sanity check: does nle() recover a posterior we already know the answer to?
##
##     Rscript dev/nle-regression-check.R [n_obs] [n_sims] [seed] [n_components] [estimator]
##
## Defaults to 500 observations, 10000 simulations, seed 1, 5 components, mdn.
##
## `estimator = "linear_gaussian"` is torch-free and finishes in seconds, and it
## is the wrong model for this problem, on purpose: see the note next to the
## nle() call. Use it to exercise the pipeline, not to check the answer.
##
## The model is a linear regression with three covariates and Gaussian noise,
## so its likelihood is available in closed form and Stan can sample the exact
## posterior. That posterior is the reference. We then throw the likelihood
## away, hand the same data to nle() as a black-box simulator, and sample the
## surrogate two ways: with the package's own slice sampler, and by exporting
## the surrogate to Stan and running NUTS on it. Three posteriors over the same
## five parameters, one of them correct by construction.
##
## Nobody should use NLE for a Gaussian regression. The point is that it is one
## of the few settings where "the learned likelihood is right" can be measured
## rather than asserted.
##
## Needs torch (libtorch) and cmdstanr with a working CmdStan.

suppressPackageStartupMessages(library(neuralsbi))

args   <- commandArgs(trailingOnly = TRUE)
n_obs  <- if (length(args) >= 1) as.integer(args[1]) else 500L
n_sim  <- if (length(args) >= 2) as.integer(args[2]) else 10000L
seed   <- if (length(args) >= 3) as.integer(args[3]) else 1L
n_comp <- if (length(args) >= 4) as.integer(args[4]) else 5L
de_kind <- if (length(args) >= 5) args[5] else "mdn"

if (de_kind != "linear_gaussian" &&
    (!requireNamespace("torch", quietly = TRUE) || !torch::torch_is_installed())) {
  stop("torch (libtorch) is required: install.packages('torch'); torch::install_torch()",
       call. = FALSE)
}
has_stan <- requireNamespace("cmdstanr", quietly = TRUE) &&
  !is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL))
if (!has_stan) {
  stop("cmdstanr with a working CmdStan is required: cmdstanr::install_cmdstan()",
       call. = FALSE)
}

set.seed(seed)
# Spread the simulator across cores if future is installed. The workers die
# with this process, so there is no plan to reset.
if (requireNamespace("future", quietly = TRUE)) future::plan(future::multisession)

## Every check below appends its verdict here, and the script exits non-zero if
## any of them failed, so this can run unattended.
checks <- list()
check <- function(label, ok, detail = "") {
  checks[[length(checks) + 1L]] <<- list(label = label, ok = isTRUE(ok),
                                         detail = detail)
  cat(sprintf("  [%s] %-42s %s\n", if (isTRUE(ok)) "PASS" else "FAIL",
              label, detail))
  invisible(ok)
}
section <- function(...) cat("\n== ", ..., "\n", sep = "")
elapsed <- function(t0) sprintf("%.1fs", as.numeric(Sys.time() - t0, units = "secs"))

# ---- the model ------------------------------------------------------------
#
# One simulator call draws one row: a set of covariates, and the response that
# goes with them. The covariates have to travel with the response because the
# regression likelihood is conditional on them and they differ from row to row.
# What the surrogate then learns is the joint density
#
#     p(y, x | theta) = p(y | x, theta) p(x),
#
# and the second factor does not involve theta, so it shifts every
# log-likelihood by the same constant and leaves the posterior alone.

rho <- 0.6                       # correlation between x1 and x2, to give the
                                 # posterior correlated coefficients

draw_covariates <- function() {
  z <- rnorm(3)
  c(x1 = z[1], x2 = rho * z[1] + sqrt(1 - rho^2) * z[2], x3 = z[3])
}

simulator <- function(b0, b1, b2, b3, sigma) {
  x <- draw_covariates()
  y <- b0 + b1 * x[1] + b2 * x[2] + b3 * x[3] + rnorm(1, 0, sigma)
  c(y = unname(y), x)
}

# sigma is kept away from zero: a residual scale near zero sends the
# conditional density of y to a spike, and the estimator would be trained on
# rows whose log-density dwarfs everything else.
prior <- prior_uniform(
  low  = c(b0 = -3, b1 = -3, b2 = -3, b3 = -3, sigma = 0.2),
  high = c(b0 =  3, b1 =  3, b2 =  3, b3 =  3, sigma = 2)
)

theta_true <- c(b0 = 0.5, b1 = 1.2, b2 = -0.8, b3 = 0.3, sigma = 0.6)
pnames <- names(theta_true)

obs <- t(replicate(n_obs, do.call(simulator, as.list(theta_true))))

section("data: ", n_obs, " rows x ", ncol(obs), " columns")
print(round(head(obs, 3), 3))

X <- cbind(1, obs[, c("x1", "x2", "x3")])
b_ols <- drop(solve(crossprod(X), crossprod(X, obs[, "y"])))
s_ols <- sqrt(sum((obs[, "y"] - X %*% b_ols)^2) / (n_obs - 4))
cat("\nOLS: ", paste(sprintf("%s=%.3f", pnames[1:4], b_ols), collapse = "  "),
    sprintf("  residual sd=%.3f\n", s_ols))

# ---- the reference: exact likelihood, in Stan -----------------------------

section("Stan, exact likelihood")
library(cmdstanr)

reg_code <- "
data {
  int<lower=1> N;
  vector[N] y;
  matrix[N, 3] X;
}
parameters {
  real<lower=-3, upper=3> b0;
  vector<lower=-3, upper=3>[3] b;
  real<lower=0.2, upper=2> sigma;
}
model {
  y ~ normal(b0 + X * b, sigma);   // flat priors, implicit in the bounds
}
"

t0 <- Sys.time()
reg_model <- cmdstan_model(write_stan_file(reg_code))
exact <- reg_model$sample(
  data = list(N = nrow(obs), y = obs[, "y"], X = obs[, c("x1", "x2", "x3")]),
  chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 2000, seed = seed, refresh = 0
)
cat("  sampled in ", elapsed(t0), "\n", sep = "")
d_exact <- as.matrix(exact$draws(c("b0", "b", "sigma"), format = "draws_matrix"))
colnames(d_exact) <- pnames
print(round(rbind(mean = colMeans(d_exact), sd = apply(d_exact, 2, sd)), 4))

# ---- the surrogate --------------------------------------------------------
#
# nle() sees the prior and the simulator. It never sees the normal density,
# the design matrix, or the fact that this is a regression.

section("nle(), ", format(n_sim, big.mark = ","), " simulations, ", de_kind)
t0 <- Sys.time()
# An MDN maps theta to a mixture over the observation and never looks at the
# observation itself, so all n rows are scored from one forward pass; a flow
# would run once per row. With Gaussian covariates and Gaussian noise the joint
# density of a row is a 4-d Gaussian, so the components are more than enough.
#
# "linear_gaussian" cannot do this problem, and it is worth knowing why. It
# fits a mean linear in the conditioning variable and one covariance shared by
# every conditioning value. Here the conditioning variable is theta, and the
# true row density given theta is N((b0, 0, 0, 0), S(theta)) with
#
#     S(theta) = [[b' Sx b + sigma^2, (Sx b)'], [Sx b, Sx]].
#
# The mean is linear in theta, so the estimator gets b0. Everything else about
# theta lives in the covariance, which the estimator holds fixed, so b1, b2, b3
# and sigma leave almost no trace in the learned likelihood and their
# posteriors fall back to the prior. Exact for a linear-Gaussian simulator
# means x = A theta + noise with the noise not depending on theta; a regression
# read as a density over (y, x) is not that model.
fit <- nle(prior, simulator, n_simulations = n_sim,
           density_estimator = de_kind, n_components = n_comp, seed = seed)
cat("  trained in ", elapsed(t0), "\n", sep = "")
print(fit)

# ---- is the learned likelihood the right function of theta? ---------------
#
# The strongest check here does not involve MCMC at all: compare the surrogate
# log-likelihood of the whole dataset against the exact one, along a slice
# through parameter space. Both are functions of theta; they should differ by a
# constant.

section("log-likelihood profile, learned vs exact")
loglik <- likelihood_fn(fit, obs)
exact_loglik <- function(theta) {
  mu <- theta[1] + obs[, c("x1", "x2", "x3")] %*% theta[2:4]
  sum(dnorm(obs[, "y"], mu, theta[5], log = TRUE))
}

grid <- seq(theta_true[["b1"]] - 0.2, theta_true[["b1"]] + 0.2, length.out = 41)
profile <- t(sapply(grid, function(b1) {
  theta <- replace(theta_true, "b1", b1)
  c(learned = loglik(theta), exact = exact_loglik(theta))
}))
profile <- sweep(profile, 2, apply(profile, 2, max))   # centred: the constant drops out

peak <- grid[apply(profile, 2, which.max)]
rms <- sqrt(mean((profile[, "learned"] - profile[, "exact"])^2))
sd_b1 <- sd(d_exact[, "b1"])
cat(sprintf("  peak at b1: exact %.3f, learned %.3f (%.2f exact posterior sd)\n",
            peak[2], peak[1], (peak[1] - peak[2]) / sd_b1))
cat(sprintf("  rms gap between the centred curves: %.2f nats\n", rms))
check("profile peaks agree within 0.5 posterior sd",
      abs(peak[1] - peak[2]) < 0.5 * sd_b1,
      sprintf("%.2f sd", (peak[1] - peak[2]) / sd_b1))

# ---- sampling it: the built-in slice sampler ------------------------------

section("nle() posterior, slice sampler")
t0 <- Sys.time()
post_slice <- posterior(fit, obs, n_chains = 20, warmup = 500, seed = seed + 1L)
d_slice <- sample(post_slice, 8000)
cat("  sampled in ", elapsed(t0), "\n", sep = "")
diag_slice <- attr(d_slice, "diagnostics")
print(diag_slice)

# ---- sampling it: the same likelihood, transpiled to Stan -----------------
#
# stan_code() writes the trained network out as Stan source with the weights
# passed as data, so NUTS differentiates the surrogate itself and nothing links
# against torch. Two samplers with nothing in common landing in the same place
# is the check that the transpiled code means what log_lik() means.

section("nle() posterior, Stan NUTS on the exported likelihood")
t0 <- Sys.time()
post_stan <- posterior(fit, obs, sampler = "stan", seed = seed + 2L,
                       iter_warmup = 1000, iter_sampling = 2000)
d_nle_stan <- sample(post_stan, 8000)
cat("  sampled in ", elapsed(t0), "\n", sep = "")
diag_stan <- attr(d_nle_stan, "diagnostics")
print(diag_stan)

check("slice sampler converged (max rhat < 1.05)",
      max(diag_slice$rhat, na.rm = TRUE) < 1.05,
      sprintf("max rhat %.3f", max(diag_slice$rhat, na.rm = TRUE)))
check("NUTS on the surrogate converged (max rhat < 1.05)",
      max(diag_stan$rhat, na.rm = TRUE) < 1.05,
      sprintf("max rhat %.3f", max(diag_stan$rhat, na.rm = TRUE)))

# ---- the comparison -------------------------------------------------------

draws <- list(exact = d_exact, slice = d_slice, `nle-stan` = d_nle_stan)

section("posterior means")
print(round(rbind(truth = theta_true, t(sapply(draws, colMeans))), 4))
section("posterior sds")
# The prior sd is the "learned nothing" reference: a posterior sitting on it is
# one the surrogate carries no information about.
print(round(rbind(prior = apply(sample_prior(prior, 20000), 2, sd),
                  t(sapply(draws, function(d) apply(d, 2, sd)))), 4))

# The units that matter are the reference posterior's own standard deviations:
# a shift of 0.1 sd is invisible in practice, a shift of 2 sd is a different
# answer quoted with the same confidence.
sd_exact <- apply(d_exact, 2, sd)
z <- t(sapply(draws[-1], function(d) (colMeans(d) - colMeans(d_exact)) / sd_exact))
section("mean shift from the exact posterior, in exact posterior sds")
print(round(z, 3))

ratio <- t(sapply(draws[-1], function(d) apply(d, 2, sd) / sd_exact))
section("posterior sd, as a ratio to the exact posterior")
print(round(ratio, 3))

section("classifier two-sample test (0.5 = indistinguishable)")
acc <- c(exact_vs_slice = c2st(d_exact, d_slice, seed = seed)$accuracy,
         exact_vs_stan  = c2st(d_exact, d_nle_stan, seed = seed)$accuracy,
         slice_vs_stan  = c2st(d_slice, d_nle_stan, seed = seed)$accuracy)
print(round(acc, 3))

# ---- how the disagreement scales with the number of rows -------------------
#
# This is what tells a broken pipeline apart from an imperfect surrogate. The
# surrogate is off by some amount per row. Summed over n rows that error grows
# like n, while the information the data carry about theta also grows like n,
# so the posterior mean settles at a fixed offset in absolute terms while the
# posterior sd shrinks like 1/sqrt(n). The shift measured in posterior sds
# therefore grows like sqrt(n), which is why `per_sqrt_n` below should be
# roughly flat down the column. A pipeline bug does not have that signature.

sweep_n <- unique(as.integer(Filter(function(k) k <= n_obs, c(25, 100, n_obs))))
if (!identical(Sys.getenv("NLE_CHECK_SWEEP"), "0") && length(sweep_n) > 1L) {
  section("agreement against exact, by number of observed rows")
  rows <- lapply(sweep_n, function(n) {
    sub <- obs[seq_len(n), , drop = FALSE]
    ex <- reg_model$sample(
      data = list(N = n, y = sub[, "y"], X = sub[, c("x1", "x2", "x3")]),
      chains = 4, parallel_chains = 4,
      iter_warmup = 1000, iter_sampling = 1000, seed = seed, refresh = 0
    )
    de <- as.matrix(ex$draws(c("b0", "b", "sigma"), format = "draws_matrix"))
    dn <- sample(posterior(fit, sub, n_chains = 20, warmup = 500,
                           seed = seed + 1L), 4000)
    zz <- (colMeans(dn) - colMeans(de)) / apply(de, 2, sd)
    c(n = n, max_abs_z = max(abs(zz)),
      per_sqrt_n = max(abs(zz)) / sqrt(n),
      c2st = c2st(de, dn, seed = seed)$accuracy)
  })
  sweep_tab <- do.call(rbind, rows)
  print(round(sweep_tab, 3))
  check("smallest dataset agrees with exact (max |shift| < 1 sd)",
        sweep_tab[1, "max_abs_z"] < 1,
        sprintf("n=%d, max |shift| %.2f sd, c2st %.2f", sweep_tab[1, "n"],
                sweep_tab[1, "max_abs_z"], sweep_tab[1, "c2st"]))
}

section("verdicts")
# The two NLE routes target the same distribution, so they should agree much
# more tightly than either agrees with the exact posterior. Measure that gap in
# the surrogate posterior's own sds, not the exact posterior's: when the
# surrogate is far wider than the reference, a gap that is negligible for the
# two samplers looks enormous in the reference's units.
gap <- max(abs(colMeans(d_slice) - colMeans(d_nle_stan)) / apply(d_slice, 2, sd))
check("slice and NUTS agree on the surrogate posterior",
      gap < 0.15 && acc[["slice_vs_stan"]] < 0.6,
      sprintf("max gap %.3f of its own sd, c2st %.2f", gap, acc[["slice_vs_stan"]]))
check("NLE posterior means within 0.5 sd of exact",
      max(abs(z)) < 0.5, sprintf("max |shift| %.2f sd", max(abs(z))))
check("NLE posterior widths within 25% of exact",
      max(abs(log(ratio))) < log(1.25),
      sprintf("sd ratio in [%.2f, %.2f]", min(ratio), max(ratio)))
check("a classifier cannot separate NLE from exact (c2st < 0.65)",
      max(acc[1:2]) < 0.65, sprintf("max c2st %.2f", max(acc[1:2])))

# ---- optional figure ------------------------------------------------------

plot_file <- Sys.getenv("NLE_CHECK_PLOT", "")
if (nzchar(plot_file)) {
  grDevices::png(plot_file, width = 1200, height = 800, res = 120)
  cols <- c(exact = "black", slice = "firebrick", `nle-stan` = "steelblue")
  op <- par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  for (j in seq_along(pnames)) {
    dens <- lapply(draws, function(d) density(d[, j]))
    plot(dens$exact, main = pnames[j], xlab = "", col = cols[1], lwd = 2,
         ylim = c(0, 1.05 * max(sapply(dens, function(d) max(d$y)))))
    for (i in 2:3) lines(dens[[i]], col = cols[i], lwd = 2, lty = i)
    abline(v = theta_true[j], col = "grey70")
  }
  plot.new()
  legend("center", names(cols), col = cols, lty = 1:3, lwd = 2, bty = "n")
  par(op)
  grDevices::dev.off()
  cat("\nwrote ", plot_file, "\n", sep = "")
}

# ---- summary --------------------------------------------------------------

failed <- Filter(function(c) !c$ok, checks)
section(sprintf("%d of %d checks passed", length(checks) - length(failed),
                length(checks)))
if (length(failed)) {
  for (f in failed) cat("  FAILED: ", f$label, " (", f$detail, ")\n", sep = "")
  cat("\nA failure here is usually the surrogate, not the sampler: the error in\n",
      "q(x | theta) is summed over ", n_obs, " rows while the information in the data\n",
      "grows at the same rate, so a small per-row error lands as a confident\n",
      "posterior in the wrong place. Raise n_simulations before anything else.\n",
      sep = "")
  quit(status = 1L)
}
cat("\nnle() reproduces the exact posterior on a problem where the exact\n",
    "posterior is known, and both samplers agree on the surrogate.\n", sep = "")
