## Does the learned likelihood converge on the right answer as the budget grows?
##
##     Rscript dev/nle-regression-budget.R [n_obs] [seed] [budgets...]
##
## Defaults to 500 observations, seed 1, budgets 2000 5000 10000 20000 40000.
##
## Companion to dev/nle-regression-check.R, which asks whether one nle() fit
## reproduces the exact posterior. This one asks the cheaper and more basic
## question: does the *maximum* of the learned likelihood walk toward the OLS
## estimate as the simulation budget grows? Estimation error shrinks with
## budget. A bug does not.
##
## There is no MCMC and no Stan here, which is what makes it quick: the
## surrogate log-likelihood is a plain function of theta, so optim() finds its
## maximum directly, and OLS is the exact maximum of the real likelihood.

suppressPackageStartupMessages(library(neuralsbi))

args    <- commandArgs(trailingOnly = TRUE)
n_obs   <- if (length(args) >= 1) as.integer(args[1]) else 500L
seed    <- if (length(args) >= 2) as.integer(args[2]) else 1L
budgets <- if (length(args) >= 3) as.integer(args[-(1:2)]) else
  c(2000L, 5000L, 10000L, 20000L, 40000L)

set.seed(seed)
if (requireNamespace("future", quietly = TRUE)) future::plan(future::multisession)

rho <- 0.6
draw_covariates <- function() {
  z <- rnorm(3)
  c(x1 = z[1], x2 = rho * z[1] + sqrt(1 - rho^2) * z[2], x3 = z[3])
}
simulator <- function(b0, b1, b2, b3, sigma) {
  x <- draw_covariates()
  y <- b0 + b1 * x[1] + b2 * x[2] + b3 * x[3] + rnorm(1, 0, sigma)
  c(y = unname(y), x)
}
prior <- prior_uniform(
  low  = c(b0 = -3, b1 = -3, b2 = -3, b3 = -3, sigma = 0.2),
  high = c(b0 =  3, b1 =  3, b2 =  3, b3 =  3, sigma = 2)
)
theta_true <- c(b0 = 0.5, b1 = 1.2, b2 = -0.8, b3 = 0.3, sigma = 0.6)
pnames <- names(theta_true)

obs <- t(replicate(n_obs, do.call(simulator, as.list(theta_true))))

X <- cbind(1, obs[, c("x1", "x2", "x3")])
b_ols <- drop(solve(crossprod(X), crossprod(X, obs[, "y"])))
s_ols <- sqrt(sum((obs[, "y"] - X %*% b_ols)^2) / n_obs)  # MLE scale, not the
                                                          # unbiased one, since
                                                          # this is an argmax
ols <- setNames(c(b_ols, s_ols), pnames)
prior_sd <- apply(sample_prior(prior, 20000), 2, sd)

cat("\nobservations : ", n_obs, "\n", sep = "")
cat("OLS / MLE    : ", paste(sprintf("%s=%.4f", pnames, ols), collapse = "  "),
    "\n", sep = "")

## The maximum of a learned likelihood, found from several starting points so a
## local optimum does not get reported as the answer. Bounds are the prior's,
## because outside them the surrogate was never trained.
argmax_loglik <- function(loglik, n_starts = 4L) {
  starts <- rbind(ols, theta_true, sample_prior(prior, n_starts))
  best <- NULL
  for (i in seq_len(nrow(starts))) {
    o <- tryCatch(
      stats::optim(starts[i, ], function(th) -loglik(th), method = "L-BFGS-B",
                   lower = prior$lower, upper = prior$upper,
                   control = list(maxit = 500)),
      error = function(e) NULL)
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
  }
  setNames(best$par, pnames)
}

rows <- lapply(budgets, function(n_sim) {
  t0 <- Sys.time()
  fit <- nle(prior, simulator, n_simulations = n_sim,
             density_estimator = "mdn", n_components = 5, seed = seed)
  secs <- as.numeric(Sys.time() - t0, units = "secs")
  est <- argmax_loglik(likelihood_fn(fit, obs))
  cat(sprintf("\n%6d simulations (%.0fs): %s\n", n_sim, secs,
              paste(sprintf("%s=%.4f", pnames, est), collapse = "  ")))
  c(n_sim = n_sim, est,
    max_abs_err = max(abs(est - ols)),
    max_rel_prior = max(abs(est - ols) / prior_sd))
})

tab <- do.call(rbind, rows)
cat("\n== maximum of the learned likelihood, against OLS on the same data\n")
print(round(tab[, c("n_sim", pnames)], 4))
cat("\n== error against OLS\n")
print(round(tab[, c("n_sim", "max_abs_err", "max_rel_prior")], 4))

# Fit a slope on the log-log scale. Estimation error falling off like a power of
# the budget is the expected shape; a flat or rising line is not.
if (nrow(tab) >= 3) {
  slope <- coef(lm(log(tab[, "max_abs_err"]) ~ log(tab[, "n_sim"])))[[2]]
  cat(sprintf("\nerror ~ budget^%.2f\n", slope))
  cat(if (slope < -0.15)
    "The learned likelihood converges on the exact one as the budget grows,\nwhich is what estimation error looks like.\n"
  else
    "The error is not falling with budget. That is not estimation error;\nlook for a bug.\n")
}
