# 07_npe_diagnostics.R -------------------------------------------------------
#
# A fitted posterior always returns draws. Nothing in the return value tells
# you whether those draws are the posterior. This script is the tour of the
# checks that do, on the fit from 01_basic_npe_example.R, and the question each
# one answers:
#
#   summary()                is the fit trained, and on how much?
#   posterior_predictive()   does the model reproduce the data it saw?
#   sbc()                    are the marginal posteriors calibrated?
#   expected_coverage()      does a 90% interval contain the truth 90% of the time?
#   tarp()                   is the JOINT posterior calibrated, not just the margins?
#   c2st()                   are two sets of draws distinguishable?
#
# 09_sbc.R, 10_c2st.R and 11_tarp.R each take one of these further.
#
# Model and data source
#   Grinsztajn, L., Semenova, E., Margossian, C. C. and Riou, J. "Bayesian
#   workflow for disease transmission modeling in Stan". Stan case study.
#   https://mc-stan.org/learn-stan/case-studies/boarding_school_case_study.html
#   The 1978 influenza A (H1N1) outbreak in a British boarding school,
#   N = 763, 14 days of boys confined to bed
#   (outbreaks::influenza_england_1978_school), fit with an SIR ODE and a
#   negative-binomial observation model. See 01_basic_npe_example.R.
#
# Runtime: about 5 minutes on a laptop CPU.

library(neuralsbi)

has_torch <- requireNamespace("torch", quietly = TRUE) &&
  torch::torch_is_installed()
has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

# ---------------------------------------------------------------------------
# The model, as in 01
# ---------------------------------------------------------------------------

cases <- c(3, 8, 26, 76, 225, 298, 258, 233, 189, 128, 68, 29, 14, 4)
n_days <- length(cases)
N <- 763
y0 <- c(S = N - 1, I = 1, R = 0)

sir_infected <- function(beta, gamma, dt = 0.25) {
  deriv <- function(y) {
    inf <- beta * y[2L] * y[1L] / N
    rec <- gamma * y[2L]
    c(-inf, inf - rec, rec)
  }
  y <- y0
  out <- numeric(n_days)
  for (s in seq_len(round(n_days / dt))) {
    k1 <- deriv(y); k2 <- deriv(y + dt / 2 * k1)
    k3 <- deriv(y + dt / 2 * k2); k4 <- deriv(y + dt * k3)
    y <- y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    t_now <- s * dt
    if (abs(t_now - round(t_now)) < 1e-8 && round(t_now) >= 1) {
      out[round(t_now)] <- max(y[2L], 0)
    }
  }
  out
}

simulator <- function(beta, gamma, phi_inv) {
  mu <- sir_infected(beta, gamma)
  stats::setNames(rnbinom(n_days, size = 1 / phi_inv, mu = mu + 1e-8),
                  paste0("day", seq_len(n_days)))
}

rtnorm0 <- function(n, m, s) {
  stats::qnorm(stats::runif(n, stats::pnorm(0, m, s), 1), m, s)
}
dtnorm0 <- function(x, m, s) {
  stats::dnorm(x, m, s, log = TRUE) -
    log(stats::pnorm(0, m, s, lower.tail = FALSE))
}

prior <- prior_custom(
  sample_fn = function(n) cbind(rtnorm0(n, 2, 1), rtnorm0(n, 0.4, 0.5),
                                stats::rexp(n, 5)),
  log_prob_fn = function(theta) {
    dtnorm0(theta[, 1], 2, 1) + dtnorm0(theta[, 2], 0.4, 0.5) +
      stats::dexp(theta[, 3], 5, log = TRUE)
  },
  dim = 3, lower = 0,
  param_names = c("beta", "gamma", "phi_inv")
)

estimator <- if (has_torch) "maf" else "linear_gaussian"
fit <- npe(prior, simulator, n_simulations = 4000,
           density_estimator = estimator, seed = 2024)

# ---------------------------------------------------------------------------
# 1. summary(): what actually got trained
# ---------------------------------------------------------------------------
#
# The cheapest check, and the one people skip. epochs_trained tells you whether
# early stopping fired or the cap did; best_val_loss is the number to watch
# when comparing architectures or simulation budgets.

info <- summary(fit)
str(info)

# ---------------------------------------------------------------------------
# 2. posterior_predictive(): does the model reproduce the data?
# ---------------------------------------------------------------------------
#
# Draw from the posterior, push those parameters back through the simulator,
# compare with the observation. This catches a fit that is confidently wrong,
# and it also catches a model that cannot produce data like yours no matter
# what the parameters are, which is a different and more serious problem.

post <- posterior(fit, x_obs = cases)
pred <- posterior_predictive(post, simulator, n = 500)

band <- apply(pred, 2, quantile, probs = c(0.05, 0.95))
cat("\ndays inside the 90% predictive band:",
    sprintf("%.0f%%", 100 * mean(cases >= band[1, ] & cases <= band[2, ])),
    "\n")

if (has_ggplot) {
  # Returns, invisibly, the fraction of predictive draws below the observation
  # in each dimension. Values near 0 or 1 mark the dimensions that misfit.
  q <- plot_posterior_predictive(pred, cases)
  print(round(q, 2))
}

# ---------------------------------------------------------------------------
# 3. sbc(): are the marginals calibrated?
# ---------------------------------------------------------------------------
#
# Draw a true parameter from the prior, simulate data from it, and rank that
# true value among posterior draws conditioned on the simulated data. Repeat.
# If the posterior is calibrated the ranks are uniform. A U-shaped histogram
# means the posterior is too narrow; an arch means it is too wide; a slope
# means it is biased.
#
# SBC needs no ground-truth posterior, which is why it works on models where
# nothing analytic exists. It tests the fit against the prior it was trained
# on, so leave `prior` at its default.

sbc_res <- sbc(fit, simulator, n_sbc = 200L, n_posterior_samples = 400L,
               seed = 1)
print(sbc_res)

# Large p-values mean the ranks look uniform. They are a screen, not a verdict.
# The test is not powerful at a few hundred trials, and it is unstable there
# too: on this fit the phi_inv p-value swings between under 0.01 and over 0.5
# depending only on n_sbc and the seed. Do not report it as the headline. The
# coverage table below is the same information on a scale that does not jump
# around, and it is what to read.

if (has_ggplot) {
  plot_sbc(sbc_res, param = "beta")
  plot_sbc(sbc_res, param = "gamma")
}

# ---------------------------------------------------------------------------
# 4. expected_coverage(): the same ranks, read as intervals
# ---------------------------------------------------------------------------
#
# A rank histogram is hard to eyeball. Coverage says the same thing in units
# people already use: of the trials, what fraction had the truth inside the
# nominal 50%, 80%, 90% central interval?

cov <- expected_coverage(sbc_res, levels = c(0.5, 0.8, 0.9, 0.95))
print(round(cov, 3))

# Expect beta and gamma to track the nominal level closely and phi_inv to sit
# a few points below it at every level. That is the useful outcome: the two
# parameters the epidemic curve identifies are calibrated, and the
# negative-binomial dispersion, which it identifies only weakly, is mildly
# overconfident. A fit can be right about what you care about and off on a
# nuisance parameter, and this table is how you tell which is which.

if (has_ggplot) plot_coverage(sbc_res)

# Above the diagonal, the posterior is too wide (conservative). Below it, too
# narrow, which is the failure that matters: it means reported intervals are
# smaller than they should be.

# ---------------------------------------------------------------------------
# 5. tarp(): is the joint calibrated?
# ---------------------------------------------------------------------------
#
# SBC ranks each parameter on its own. A posterior can have perfect marginals
# and the wrong correlation between beta and gamma, and SBC will not see it.
# TARP measures distances in the full parameter space against random reference
# points, so it does.

tarp_res <- tarp(fit, simulator, n_tarp = 200L, n_posterior_samples = 400L,
                 seed = 1)
print(tarp_res)

if (has_ggplot) plot_tarp(tarp_res)

# ---------------------------------------------------------------------------
# 6. c2st(): are two sets of draws distinguishable?
# ---------------------------------------------------------------------------
#
# This model has no analytic posterior to compare against, so use c2st() the
# other way: as an agreement check between two fits that should agree. If two
# independently trained estimators disagree, at least one of them is wrong.

fit_b <- npe(prior, simulator, n_simulations = 4000,
             density_estimator = estimator, patience = 12L, seed = 99)
draws_a <- sample(post, 3000)
draws_b <- sample(posterior(fit_b, x_obs = cases), 3000)
print(c2st(draws_a, draws_b, seed = 1))

# The classifier here is linear, so it sees a shift in location easily and is
# close to blind to two sets that share a mean and differ in spread. Read the
# number alongside the moments, not instead of them.
print(round(rbind(a = colMeans(draws_a), b = colMeans(draws_b)), 3))
print(round(rbind(a = apply(draws_a, 2, sd), b = apply(draws_b, 2, sd)), 3))

# Expect these two to be distinguishable. beta and gamma agree closely; the
# gap is in phi_inv, the same parameter SBC flagged, and it shows up as
# different spreads rather than different centres. Two runs disagreeing on a
# parameter's width at a fixed budget is the honest reading of how much that
# width is worth: report it, and raise n_simulations if it matters.

# ---------------------------------------------------------------------------
# 7. What a failing diagnostic usually means
# ---------------------------------------------------------------------------
#
# Ranks pile up at the edges (posterior too narrow):
#   the usual cause is too small a simulation budget. Raise n_simulations
#   before touching the architecture.
#
# Ranks pile up in the middle (posterior too wide):
#   often an underfit network. Check epochs_trained; if early stopping fired at
#   the patience limit with the loss still falling, raise patience.
#
# One parameter is calibrated and another is not:
#   look at whether the data identify it at all. A parameter the simulator
#   barely responds to will have a posterior close to the prior, and that is
#   the model being honest, not the estimator failing.
#
# TARP fails while SBC passes:
#   the marginals are right and the dependence is not. A more flexible
#   estimator ("nsf" instead of "mdn") or more simulations is the first thing
#   to try.
