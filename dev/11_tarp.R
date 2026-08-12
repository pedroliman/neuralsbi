# 11_tarp.R ------------------------------------------------------------------
#
# TARP: tests of accuracy with random points. SBC ranks each parameter on its
# own, so a posterior with perfect marginals and the wrong dependence between
# parameters passes it. TARP measures distances in the full parameter space
# against random reference points, so it is a joint test and does not have that
# blind spot.
#
# Method
#   Lemos, P., Coogan, A., Hezaveh, Y. and Perreault-Levasseur, L. (2023),
#   "Sampling-based accuracy testing of posterior estimators for general
#   inference", ICML. arXiv:2302.03026.
#
#   For each trial: draw theta from the prior, simulate x, sample the posterior
#   at x, and pick a random reference point. The fraction of posterior draws
#   closer to the reference than the truth is the credibility level of the
#   smallest distance-based credible region containing the truth. For a
#   calibrated posterior those fractions are uniform, so the expected coverage
#   probability at level alpha equals alpha.
#
# Task source
#   Lueckmann, J.-M. et al. "Benchmarking Simulation-Based Inference",
#   AISTATS 2021. https://github.com/sbi-benchmark/sbibm
#   task_two_moons() is sbibm's two_moons: uniform prior on [-1, 1]^2, and a
#   simulator that maps parameters onto a noisy crescent. The posterior for a
#   typical observation has two symmetric modes joined by a curved ridge, which
#   is exactly the structure a marginal test cannot see.
#
# Runtime: about 3 minutes on a laptop CPU.

library(neuralsbi)

has_torch <- requireNamespace("torch", quietly = TRUE) &&
  torch::torch_is_installed()
has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

task <- task_two_moons()
print(task)
prior <- task$prior
simulator <- task$simulator

# ---------------------------------------------------------------------------
# A fit that should do well, and one that cannot
# ---------------------------------------------------------------------------
#
# A flow can represent a curved bimodal posterior. The closed-form conditional
# Gaussian cannot: it has one mode and elliptical contours, by construction.
# Running both makes the difference between the two diagnostics visible.

fit_flow <- if (has_torch) {
  npe(prior, simulator, n_simulations = 5000, density_estimator = "nsf",
      seed = 1)
} else NULL

fit_gauss <- npe(prior, simulator, n_simulations = 5000,
                 density_estimator = "linear_gaussian", seed = 1)

# ---------------------------------------------------------------------------
# What the two posteriors look like
# ---------------------------------------------------------------------------

set.seed(4)
theta_true <- sample_prior(prior, 1)
x_obs <- simulator(theta_true[1, ])

d_gauss <- sample(posterior(fit_gauss, x_obs = x_obs), 3000)
cat("linear_gaussian posterior:\n"); print(summary(d_gauss))

if (!is.null(fit_flow)) {
  d_flow <- sample(posterior(fit_flow, x_obs = x_obs), 3000)
  cat("\nnsf posterior:\n"); print(summary(d_flow))
  cat("\ncorrelation between the two parameters:\n")
  cat(sprintf("  linear_gaussian %.3f   nsf %.3f\n",
              cor(d_gauss[, 1], d_gauss[, 2]),
              cor(d_flow[, 1], d_flow[, 2])))
  if (has_ggplot && requireNamespace("GGally", quietly = TRUE) &&
      requireNamespace("ggdensity", quietly = TRUE)) {
    pairplot(d_flow, truth = theta_true[1, ])
    pairplot(d_gauss, truth = theta_true[1, ])
  }
}

# ---------------------------------------------------------------------------
# TARP
# ---------------------------------------------------------------------------

tarp_gauss <- tarp(fit_gauss, simulator, n_tarp = 300L,
                   n_posterior_samples = 500L, seed = 1)
print(tarp_gauss)

if (!is.null(fit_flow)) {
  tarp_flow <- tarp(fit_flow, simulator, n_tarp = 300L,
                    n_posterior_samples = 500L, seed = 1)
  print(tarp_flow)
}

# print() reports max |ECP - nominal|, which is the single number to compare
# fits by. Zero is perfect calibration.

# ---------------------------------------------------------------------------
# Reading the curve
# ---------------------------------------------------------------------------

str(tarp_gauss[c("levels", "ecp", "n_tarp", "references")])

if (has_ggplot) {
  plot_tarp(tarp_gauss)
  if (!is.null(fit_flow)) plot_tarp(tarp_flow)
}

# Above the diagonal: the posterior is too wide, and reported regions contain
# the truth more often than advertised. Below: too narrow, and they contain it
# less often, which is the failure that costs you. The shaded band is the
# Monte-Carlo uncertainty from a finite n_tarp, so a curve inside the band is
# not evidence of anything.

# ---------------------------------------------------------------------------
# Reference points
# ---------------------------------------------------------------------------
#
#   "uniform" (default) draws references uniformly over the hyper-rectangle
#     spanned by the true parameter draws, as in the paper.
#   "prior" draws them from the prior.
#
# They differ when the prior is far from uniform: prior references concentrate
# where the prior does, so the test looks hardest where most of the prior mass
# is. On this task the prior is uniform on a box, so the two should agree
# closely, which is a useful way to see that the choice is not doing the work.

tarp_prior_ref <- tarp(fit_gauss, simulator, n_tarp = 300L,
                       n_posterior_samples = 500L, references = "prior",
                       seed = 1)
print(data.frame(
  nominal = tarp_gauss$levels,
  uniform_refs = round(tarp_gauss$ecp, 3),
  prior_refs = round(tarp_prior_ref$ecp, 3)
))

# ---------------------------------------------------------------------------
# TARP and SBC together
# ---------------------------------------------------------------------------
#
# Run both. They cost about the same and they fail on different things.
# Distances in TARP are computed after z-scoring each parameter by the spread
# of the true draws, so parameters on different scales contribute comparably;
# nothing similar is needed in SBC because it never mixes parameters.

sbc_gauss <- sbc(fit_gauss, simulator, n_sbc = 300L,
                 n_posterior_samples = 500L, seed = 1)
print(sbc_gauss)
cat(sprintf("TARP max deviation, linear_gaussian: %.3f\n",
            max(abs(tarp_gauss$ecp - tarp_gauss$levels))))

if (!is.null(fit_flow)) {
  sbc_flow <- sbc(fit_flow, simulator, n_sbc = 300L,
                  n_posterior_samples = 500L, seed = 1)
  print(sbc_flow)
  cat(sprintf("TARP max deviation, nsf: %.3f\n",
              max(abs(tarp_flow$ecp - tarp_flow$levels))))
}

# On a bimodal, curved posterior the elliptical estimator has to cover both
# modes with one blob. That inflates the regions, so expect it to look
# conservative rather than overconfident: ECP above the diagonal, and marginal
# ranks that arch rather than form a U.

# ---------------------------------------------------------------------------
# Cost
# ---------------------------------------------------------------------------
#
# n_tarp simulations and n_tarp posterior samplings, the same bill as sbc().
# Cheap on an NPE fit, expensive on NLE or NRE where every trial is its own
# MCMC run; the MCMC controls reach posterior() through ...:
#
#   tarp(nle_fit, simulator, n_tarp = 40L, n_chains = 10, warmup = 100)
#
# A trial whose simulation returns non-finite output is dropped and lowers the
# effective n_tarp, which print() reports. A trial whose posterior comes back
# short is an error, for the same reason as in sbc(): it would not be
# comparable with the rest.
