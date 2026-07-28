## Why is the surrogate learning the coefficients so slowly?
##
##     Rscript dev/nle-regression-conditional.R [n_obs] [n_simulations] [seed]
##
## Defaults to 500 observations, 5000 simulations, seed 1.
##
## The suspicion this script tests: it is not that the covariates fail to vary
## (they are redrawn on every simulator call), it is the *role* they play.
##
## nle() learns the density of a whole observation given theta. Here the
## observation is the row (y, x1, x2, x3), so the surrogate learns the joint
##
##     p(y, x | theta) = N( (b0, 0, 0, 0),  S(theta) ),
##     S(theta) = [[b' Sx b + sigma^2, (Sx b)'], [Sx b, Sx]].
##
## Look at where each parameter went. b0 sits in the mean. b1, b2, b3 appear
## nowhere in the mean at all: their only trace is the covariance between y and
## the covariates. A mean is an easy thing for a network to learn precisely; a
## covariance function of theta is not, and the coefficients are read out of it
## by inverting Sx. On top of that the estimator spends capacity on p(x), which
## carries no information about theta whatsoever.
##
## The alternative, and the point of the comparison below, is to hand the
## covariates over as conditioning variables instead of asking the estimator to
## model them. Train on q(y | x, theta) by putting the covariates into the
## conditioning vector alongside the parameters. Then the target is scalar, and
## the coefficients are back in the mean, where b0 already is:
##
##     y | x, theta ~ N(b0 + b1 x1 + b2 x2 + b3 x3, sigma^2).
##
## Both fits get the same simulation budget. Neither uses MCMC: the maximum of
## the learned likelihood is found with optim(), and OLS is the exact maximum of
## the real likelihood, so the gap between them is the surrogate's error.
##
## Note that the conditional fit is a diagnostic, not a workflow. posterior()
## would treat those extra conditioning slots as parameters to sample, so there
## is no supported way to run MCMC on it today. What it establishes is where the
## difficulty actually lives.

suppressPackageStartupMessages(library(neuralsbi))

args  <- commandArgs(trailingOnly = TRUE)
n_obs <- if (length(args) >= 1) as.integer(args[1]) else 500L
n_sim <- if (length(args) >= 2) as.integer(args[2]) else 5000L
seed  <- if (length(args) >= 3) as.integer(args[3]) else 1L

set.seed(seed)
if (requireNamespace("future", quietly = TRUE)) future::plan(future::multisession)

rho <- 0.6
draw_covariates <- function() {
  z <- rnorm(3)
  c(x1 = z[1], x2 = rho * z[1] + sqrt(1 - rho^2) * z[2], x3 = z[3])
}
theta_true <- c(b0 = 0.5, b1 = 1.2, b2 = -0.8, b3 = 0.3, sigma = 0.6)
pnames <- names(theta_true)

low  <- c(b0 = -3, b1 = -3, b2 = -3, b3 = -3, sigma = 0.2)
high <- c(b0 =  3, b1 =  3, b2 =  3, b3 =  3, sigma = 2)

# ---- the data, identical to the other two scripts -------------------------

simulator_joint <- function(b0, b1, b2, b3, sigma) {
  x <- draw_covariates()
  y <- b0 + b1 * x[1] + b2 * x[2] + b3 * x[3] + rnorm(1, 0, sigma)
  c(y = unname(y), x)
}
prior_joint <- prior_uniform(low = low, high = high)

obs <- t(replicate(n_obs, do.call(simulator_joint, as.list(theta_true))))
X <- cbind(1, obs[, c("x1", "x2", "x3")])
b_ols <- drop(solve(crossprod(X), crossprod(X, obs[, "y"])))
ols <- setNames(c(b_ols, sqrt(sum((obs[, "y"] - X %*% b_ols)^2) / n_obs)), pnames)
prior_sd <- apply(sample_prior(prior_joint, 20000), 2, sd)

cat("\nobservations : ", n_obs, "\nbudget       : ", n_sim, "\n", sep = "")
cat("OLS / MLE    : ", paste(sprintf("%s=%.4f", pnames, ols), collapse = "  "),
    "\n", sep = "")

# ---- fit A: the covariates are part of what gets modelled -----------------

cat("\n== A: q(y, x1, x2, x3 | theta), the covariates modelled\n")
t0 <- Sys.time()
fit_joint <- nle(prior_joint, simulator_joint, n_simulations = n_sim,
                 density_estimator = "mdn", n_components = 5, seed = seed)
cat("   trained in ", sprintf("%.0fs", as.numeric(Sys.time() - t0, units = "secs")),
    "\n", sep = "")
loglik_joint <- likelihood_fn(fit_joint, obs)

# ---- fit B: the covariates are conditioned on -----------------------------
#
# The conditioning vector is (theta, x) and the target is y alone. Training
# draws come from the prior for the parameter block and from the covariate
# distribution for the rest, which is the same distribution the observed
# covariates came from. That matters: the surrogate is only trained where the
# conditioning vectors were placed.

cat("\n== B: q(y | x1, x2, x3, theta), the covariates conditioned on\n")
simulator_cond <- function(b0, b1, b2, b3, sigma, x1, x2, x3) {
  c(y = b0 + b1 * x1 + b2 * x2 + b3 * x3 + rnorm(1, 0, sigma))
}
prior_cond <- prior_custom(
  sample_fn = function(n) {
    th <- sweep(matrix(runif(n * 5), n, 5), 2, high - low, `*`)
    th <- sweep(th, 2, low, `+`)
    xs <- t(replicate(n, draw_covariates()))
    out <- cbind(th, xs)
    colnames(out) <- c(pnames, "x1", "x2", "x3")
    out
  },
  dim = 8L,
  # Bounds only fence in the parameter block; the covariates are unbounded and
  # get generous limits so within_support() never rejects a real one.
  lower = c(low, x1 = -8, x2 = -8, x3 = -8),
  upper = c(high, x1 = 8, x2 = 8, x3 = 8)
)

t0 <- Sys.time()
fit_cond <- nle(prior_cond, simulator_cond, n_simulations = n_sim,
                density_estimator = "mdn", n_components = 5, seed = seed)
cat("   trained in ", sprintf("%.0fs", as.numeric(Sys.time() - t0, units = "secs")),
    "\n", sep = "")

# Each row of the dataset is scored against its own covariates, so the
# likelihood is the diagonal of the theta-by-observation cross product rather
# than a plain i.i.d. sum.
y_obs <- obs[, "y", drop = FALSE]
xs_obs <- obs[, c("x1", "x2", "x3")]
loglik_cond <- function(theta) {
  th <- cbind(matrix(theta, nrow = n_obs, ncol = 5L, byrow = TRUE), xs_obs)
  sum(diag(log_lik(fit_cond, theta = th, x = y_obs, sum_iid = FALSE)))
}

# ---- where does each one put its maximum? ---------------------------------

argmax_loglik <- function(loglik, n_starts = 3L) {
  starts <- rbind(ols, theta_true, sample_prior(prior_joint, n_starts))
  best <- NULL
  for (i in seq_len(nrow(starts))) {
    o <- tryCatch(
      stats::optim(starts[i, ], function(th) -loglik(th), method = "L-BFGS-B",
                   lower = low, upper = high, control = list(maxit = 300)),
      error = function(e) NULL)
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
  }
  setNames(best$par, pnames)
}

t0 <- Sys.time()
est_joint <- argmax_loglik(loglik_joint)
est_cond <- argmax_loglik(loglik_cond)
cat("\n   optimized in ",
    sprintf("%.0fs", as.numeric(Sys.time() - t0, units = "secs")), "\n", sep = "")

tab <- rbind(truth = theta_true, ols = ols,
             `A: covariates modelled` = est_joint,
             `B: covariates conditioned on` = est_cond)
cat("\n== maximum of the learned likelihood\n")
print(round(tab, 4))

err <- rbind(`A: covariates modelled` = est_joint - ols,
             `B: covariates conditioned on` = est_cond - ols)
cat("\n== error against OLS\n")
print(round(err, 4))
cat("\n== error against OLS, as a fraction of the prior sd\n")
print(round(sweep(abs(err), 2, prior_sd, "/"), 4))

# b0 enters the mean of the joint density; b1, b2, b3 enter only its
# covariance. If that distinction is what makes the coefficients hard, fit A
# should be much worse on the slopes than on the intercept, and fit B should
# not care.
cat("\n== intercept versus slopes (b0 is in the mean either way;\n",
    "   b1, b2, b3 are in the covariance for A and in the mean for B)\n", sep = "")
split <- cbind(b0 = abs(err[, "b0"]),
               slopes = apply(abs(err[, c("b1", "b2", "b3")]), 1, max),
               ratio = apply(abs(err[, c("b1", "b2", "b3")]), 1, max) /
                 abs(err[, "b0"]))
print(round(split, 4))

cat("\n")
if (max(abs(err[2, ])) < max(abs(err[1, ]))) {
  cat(sprintf(paste0("Conditioning on the covariates cuts the worst error from",
                     " %.3f to %.3f at the same\nbudget of %d simulations.",
                     " The budget is not the whole story: what the\nestimator",
                     " is asked to model is.\n"),
              max(abs(err[1, ])), max(abs(err[2, ])), n_sim))
} else {
  cat(sprintf(paste0("Conditioning on the covariates did not help (worst error",
                     " %.3f against %.3f).\nThe cost here is sample complexity,",
                     " not the joint formulation.\n"),
              max(abs(err[2, ])), max(abs(err[1, ]))))
}
