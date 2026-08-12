# 04_nle.R -------------------------------------------------------------------
#
# Neural likelihood estimation. Where npe() learns p(theta | x) directly,
# nle() learns a surrogate likelihood q(x | theta) and gets the posterior from
# Bayes' rule with MCMC. The payoff: the surrogate is trained on ONE
# observation, so the log-likelihood of n independent observations is a sum of
# n evaluations. Train once, then condition on 1 subject or 12 without
# retraining.
#
# Model and data source
#   Pinheiro, J. C. and Bates, D. M. "Mixed-Effects Models in S and S-PLUS"
#   (Springer, 2000), chapter 8. The first-order open-compartment model is
#   nlme::SSfol; the fit is
#
#     nlme(conc ~ SSfol(Dose, Time, lKe, lKa, lCl), data = Theoph,
#          fixed = lKe + lKa + lCl ~ 1,
#          random = pdDiag(lKe + lKa + lCl ~ 1))
#
#   see ?nlme::Theoph and ?nlme::SSfol. The data are datasets::Theoph, shipped
#   with base R: serum theophylline concentrations in 12 subjects over 25 hours
#   after a single oral dose, from Boeckmann, Sheiner and Beal (1994), NONMEM
#   Users Guide.
#
#   Two simplifications, both stated rather than hidden. Doses differ by
#   subject (3.10 to 5.86 mg/kg), so we work with dose-normalized
#   concentrations, which is exact for this dose-proportional model and makes
#   the subjects exchangeable. Sampling times differ by a few minutes between
#   subjects; we use the across-subject mean time for each of the 10 post-dose
#   draws. Pinheiro and Bates give each random effect its own variance; we use
#   one shared between-subject SD, omega, so the parameter vector stays small.
#
# Runtime: about 3 minutes on a laptop CPU.

library(neuralsbi)

has_torch <- requireNamespace("torch", quietly = TRUE) &&
  torch::torch_is_installed()

# ---------------------------------------------------------------------------
# The data
# ---------------------------------------------------------------------------

theoph <- datasets::Theoph
by_subject <- split(theoph, theoph$Subject)

# Drop the pre-dose draw: the model predicts exactly 0 there, so it carries no
# information and would give the density estimator a degenerate dimension.
keep <- 2:11
times <- colMeans(do.call(rbind, lapply(by_subject, function(s) s$Time)))[keep]
dose  <- vapply(by_subject, function(s) s$Dose[1], numeric(1))

# 12 subjects x 10 times, dose-normalized.
x_obs <- do.call(rbind, lapply(seq_along(by_subject), function(i) {
  by_subject[[i]]$conc[keep] / dose[i]
}))
colnames(x_obs) <- paste0("t", round(times, 2))
print(round(x_obs, 3))

# ---------------------------------------------------------------------------
# The simulator
# ---------------------------------------------------------------------------
#
# One call returns one subject's dose-normalized concentration profile:
#
#   c(t) / Dose = exp(lKe + lKa - lCl) *
#                 (exp(-exp(lKe) t) - exp(-exp(lKa) t)) / (exp(lKa) - exp(lKe))
#
# with a subject-level random effect of SD omega on each of the three log
# parameters and additive residual error of SD sigma.

simulator <- function(lKe, lKa, lCl, log_omega, log_sigma) {
  omega <- exp(log_omega)
  b <- rnorm(3, sd = omega)
  ke <- exp(lKe + b[1])
  ka <- exp(lKa + b[2])
  cl <- exp(lCl + b[3])
  # ka == ke is a measure-zero coincidence, but nudge it so the denominator
  # never blows up on a draw that happens to land there.
  if (abs(ka - ke) < 1e-6) ka <- ka + 1e-6
  mu <- ke * ka / (cl * (ka - ke)) * (exp(-ke * times) - exp(-ka * times))
  stats::setNames(mu + rnorm(length(times), sd = exp(log_sigma)),
                  colnames(x_obs))
}

# Ranges around the population estimates Pinheiro and Bates report
# (lKe about -2.43, lKa about 0.45, lCl about -3.21).
prior <- prior_uniform(
  low  = c(lKe = -3.5, lKa = -1.0, lCl = -4.0,
           log_omega = log(0.05), log_sigma = log(0.01)),
  high = c(lKe = -1.5, lKa =  2.0, lCl = -2.5,
           log_omega = log(0.80), log_sigma = log(0.40))
)
print(prior)

# ---------------------------------------------------------------------------
# Fit the surrogate likelihood
# ---------------------------------------------------------------------------
#
# nle() runs the same training loop as npe() with the roles of theta and x
# swapped: the estimator's target is the data and it conditions on the
# parameters. "mdn" is a good default here because its i.i.d. evaluation path
# is specialized: the mixture depends on theta alone, so it is built once and
# scored against every observation.

estimator <- if (has_torch) "mdn" else "linear_gaussian"

fit <- nle(prior, simulator, n_simulations = 4000,
           density_estimator = estimator, seed = 2024, verbose = TRUE)
print(fit)

# Note what dim_x says: 10, the width of ONE observation, not 120.

# ---------------------------------------------------------------------------
# log_lik(): the surrogate as a likelihood
# ---------------------------------------------------------------------------

theta_pb <- c(lKe = -2.43, lKa = 0.45, lCl = -3.21,
              log_omega = log(0.2), log_sigma = log(0.05))

# All 12 subjects, summed.
cat("\nlog-lik at Pinheiro-Bates values, 12 subjects:",
    sprintf("%.1f", log_lik(fit, theta_pb, x_obs)), "\n")

# One subject only.
cat("log-lik, subject 1 alone:            ",
    sprintf("%.1f", log_lik(fit, theta_pb, x_obs[1, , drop = FALSE])), "\n")

# sum_iid = FALSE returns the per-observation terms instead of the sum.
per_subject <- log_lik(fit, theta_pb, x_obs, sum_iid = FALSE)
cat("per-subject terms:", sprintf("%.1f", per_subject), "\n")

# Vectorized over rows of theta.
grid <- cbind(lKe = seq(-3.2, -1.8, length.out = 7),
              lKa = 0.45, lCl = -3.21,
              log_omega = log(0.2), log_sigma = log(0.05))
print(data.frame(lKe = grid[, "lKe"],
                 loglik = round(log_lik(fit, grid, x_obs), 1)))

# ---------------------------------------------------------------------------
# likelihood_fn(): the surrogate as a plain R function
# ---------------------------------------------------------------------------
#
# Fix the observation and you have function(theta) that anything in R can use:
# optim(), an MCMC package of your choice, an importance sampler.

loglik <- likelihood_fn(fit, x_obs)

# prior_uniform() carries its box as $lower/$upper, which is exactly what
# L-BFGS-B wants.
opt <- stats::optim(theta_pb, function(th) -loglik(th),
                    method = "L-BFGS-B",
                    lower = prior$lower, upper = prior$upper)
cat("\nsurrogate MLE:\n")
print(round(opt$par, 3))
cat("Pinheiro-Bates fixed effects: lKe -2.43, lKa 0.45, lCl -3.21\n")

# ---------------------------------------------------------------------------
# The posterior, by MCMC
# ---------------------------------------------------------------------------
#
# An NLE posterior has no closed form and no direct sampler, so posterior()
# returns an object that samples with MCMC. The default is a vectorized slice
# sampler: nothing to tune, no dependencies, and every chain advances in one
# batched call per step, so more chains cost almost nothing.

t0 <- Sys.time()
post <- posterior(fit, x_obs, n_chains = 10, warmup = 150, thin = 2, seed = 5)
draws <- sample(post, 2000)
cat(sprintf("\nMCMC: %.1f s\n", as.numeric(Sys.time() - t0, units = "secs")))

print(post)
print(summary(draws))

# ---------------------------------------------------------------------------
# The point of NLE: conditioning on however many observations you have
# ---------------------------------------------------------------------------
#
# The same trained network, conditioned on 1, 4 and 12 subjects. An NPE fit
# would need retraining for each, because its input width is fixed at training
# time.

for (n_sub in c(1, 4, 12)) {
  p <- posterior(fit, x_obs[seq_len(n_sub), , drop = FALSE],
                 n_chains = 10, warmup = 150, thin = 2, seed = 5)
  d <- sample(p, 1000)
  cat(sprintf("%2d subject(s): lCl = %.3f (sd %.3f)\n",
              n_sub, mean(d[, "lCl"]), sd(d[, "lCl"])))
}

# The posterior should tighten roughly as 1/sqrt(n) in the parameters the data
# speak to. If it does not move at all, the surrogate is not seeing that
# parameter, which is worth knowing before you trust it.

# ---------------------------------------------------------------------------
# log_prob() on an MCMC posterior
# ---------------------------------------------------------------------------
#
# It is available, and it is unnormalized: the surrogate gives the likelihood
# and the prior gives the rest, but the evidence is not computed.

print(round(log_prob(post, rbind(theta_pb, opt$par)), 2))

# ---------------------------------------------------------------------------
# Saving
# ---------------------------------------------------------------------------

path <- tempfile(fileext = ".rds")
save_nle(fit, path)
fit2 <- load_nle(path)
cat("\nreloaded fit, same log-lik:",
    isTRUE(all.equal(log_lik(fit, theta_pb, x_obs),
                     log_lik(fit2, theta_pb, x_obs))), "\n")
unlink(path)

# ---------------------------------------------------------------------------
# When to use NLE and when not to
# ---------------------------------------------------------------------------
#
# NLE, when the observation is n exchangeable units and n varies, or when the
# learned likelihood has to become one term in a larger model you write by
# hand (see 08_nle_with_stan.R).
#
# NPE, when there is one fixed observation of fixed width and you want draws in
# a forward pass instead of an MCMC run.
#
# NRE, when the data are discrete, mixed-type or otherwise awkward to model as
# a density but easy to discriminate (see 05_nre.R).
