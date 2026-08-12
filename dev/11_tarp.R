# 11_tarp.R ------------------------------------------------------------------
#
# TARP: tests of accuracy with random points. sbc() ranks each parameter on its
# own. TARP measures distances in the full parameter space against random
# reference points, so it scores the joint posterior instead of the margins.
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
#   probability (ECP) at credibility level alpha equals alpha.
#
# Task source
#   Lueckmann, J.-M., Boelts, J., Greenberg, D., Goncalves, P. and Macke, J.
#   "Benchmarking Simulation-Based Inference", AISTATS 2021.
#   https://github.com/sbi-benchmark/sbibm
#   task_gaussian_linear() (conjugate, analytic posterior) is the case where
#   everything should pass. task_two_moons() (uniform prior on [-1, 1]^2, a
#   curved bimodal posterior) is the case where a conditional-Gaussian
#   estimator cannot be right no matter how many simulations it gets.
#
# Runtime: about 3 minutes on a laptop CPU.

library(neuralsbi)

has_torch <- requireNamespace("torch", quietly = TRUE) &&
  torch::torch_is_installed()
has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

# ---------------------------------------------------------------------------
# A fit that should pass
# ---------------------------------------------------------------------------

gl <- task_gaussian_linear(dim = 3L)
print(gl)

fit_good <- npe(gl$prior, gl$simulator, n_simulations = 4000,
                density_estimator = if (has_torch) "maf" else "linear_gaussian",
                max_epochs = 150L, patience = 10L, seed = 1)

tarp_good <- tarp(fit_good, gl$simulator, n_tarp = 300L,
                  n_posterior_samples = 400L, seed = 1)
print(tarp_good)

# print() reports max |ECP - nominal|, which is the single number to compare
# fits by. Zero is perfect calibration, and at 300 trials the Monte-Carlo floor
# is around 0.05, so anything under that is noise.

# ---------------------------------------------------------------------------
# Reading the object
# ---------------------------------------------------------------------------

str(tarp_good[c("levels", "ecp", "n_tarp", "n_posterior_samples",
                "references")])

# coverage_values holds the per-trial fractions. Under calibration they are
# uniform on (0, 1), and the ECP curve is just their empirical CDF.
print(round(quantile(tarp_good$coverage_values, seq(0, 1, 0.25)), 3))

if (has_ggplot) plot_tarp(tarp_good)

# Above the diagonal: the posterior is too wide, and reported regions contain
# the truth more often than advertised. Below: too narrow, and they contain it
# less often, which is the failure that costs you. The shaded band is the
# Monte-Carlo uncertainty from a finite n_tarp, so a curve inside the band is
# not evidence of anything.

# ---------------------------------------------------------------------------
# A fit that cannot be right
# ---------------------------------------------------------------------------
#
# two_moons has a curved, bimodal posterior. "linear_gaussian" is a
# closed-form conditional Gaussian: one mode, elliptical contours. There is no
# simulation budget that fixes that.

tm <- task_two_moons()
fit_wrong <- npe(tm$prior, tm$simulator, n_simulations = 4000,
                 density_estimator = "linear_gaussian", seed = 1)

set.seed(4)
theta_true <- sample_prior(tm$prior, 1)
x_obs <- tm$simulator(theta_true[1, ])
d_wrong <- sample(posterior(fit_wrong, x_obs = x_obs), 3000)
print(summary(d_wrong))
cat(sprintf("correlation between the two parameters: %.3f\n",
            cor(d_wrong[, 1], d_wrong[, 2])))

if (has_ggplot && requireNamespace("GGally", quietly = TRUE) &&
    requireNamespace("ggdensity", quietly = TRUE)) {
  pairplot(d_wrong, truth = theta_true[1, ])
}

tarp_wrong <- tarp(fit_wrong, tm$simulator, n_tarp = 300L,
                   n_posterior_samples = 400L, seed = 1)
print(tarp_wrong)
if (has_ggplot) plot_tarp(tarp_wrong)

# ---------------------------------------------------------------------------
# Run both tests, because neither one dominates
# ---------------------------------------------------------------------------
#
# TARP is motivated as the test that sees what marginal ranks cannot. That
# motivation is real, and it does not make TARP the sharper test in every case.
# On this fit both of them fire, and SBC makes the size of the failure much
# more obvious.

sbc_wrong <- sbc(fit_wrong, tm$simulator, n_sbc = 300L,
                 n_posterior_samples = 400L, seed = 1)
print(sbc_wrong)
print(round(expected_coverage(sbc_wrong, c(0.5, 0.8, 0.9, 0.95)), 3))
cat(sprintf("TARP max |ECP - nominal|: %.3f\n",
            max(abs(tarp_wrong$ecp - tarp_wrong$levels))))

# Expect marginal coverage far below nominal, near 0.35 where 0.50 was
# promised, against a TARP deviation of about 0.10: detected, but only twice
# the Monte-Carlo floor the calibrated fit sat at.
#
# The reason is geometric. The fitted Gaussian is a long ellipse stretched
# along the ridge that joins the two crescents, so a distance-based region
# around a random reference point is not far off the right size. Project that
# same ellipse onto either axis and it is much too narrow, which is exactly
# what a rank histogram measures.
#
# The lesson is not that TARP is weak. It is that a diagnostic sees the failure
# it is shaped to see, so run both: they cost the same and they respond to
# different things. 09_sbc.R adds a third case, a posterior that ignores the
# data entirely and passes both.

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
# is. two_moons has a uniform prior on a box, so the two should agree closely
# here, which is a useful way to see that the choice is not doing the work.

tarp_prior_ref <- tarp(fit_wrong, tm$simulator, n_tarp = 300L,
                       n_posterior_samples = 400L, references = "prior",
                       seed = 1)
print(data.frame(nominal = tarp_wrong$levels,
                 uniform_refs = round(tarp_wrong$ecp, 3),
                 prior_refs = round(tarp_prior_ref$ecp, 3)))

# On the Gaussian task the prior is N(0, 0.1 I) rather than uniform, so there
# the two choices really can differ.
print(data.frame(
  nominal = tarp_good$levels,
  uniform_refs = round(tarp_good$ecp, 3),
  prior_refs = round(tarp(fit_good, gl$simulator, n_tarp = 300L,
                          n_posterior_samples = 400L, references = "prior",
                          seed = 1)$ecp, 3)))

# ---------------------------------------------------------------------------
# Scaling
# ---------------------------------------------------------------------------
#
# Distances are computed after z-scoring each parameter by the spread of the
# true draws, so a parameter measured in thousands does not swamp one measured
# in tenths. Nothing similar is needed in sbc(), which never mixes parameters.

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
