# 09_sbc.R -------------------------------------------------------------------
#
# Simulation-based calibration. The one check that works on any model, needs no
# ground-truth posterior, and tells you something a good-looking posterior
# cannot: whether the intervals mean what they say.
#
# The idea, from Talts, Betancourt, Simpson, Vehtari and Gelman (2018),
# "Validating Bayesian Inference Algorithms with Simulation-Based Calibration",
# arXiv:1804.06788. Draw theta from the prior, simulate x from it, then rank
# theta among posterior draws conditioned on that x. Averaged over draws from
# the prior, the posterior IS the prior, so those ranks are uniform whenever
# the posterior is correct. Any departure from uniformity is a real defect.
#
# Task source
#   Lueckmann, J.-M., Boelts, J., Greenberg, D., Goncalves, P. and Macke, J.
#   "Benchmarking Simulation-Based Inference", AISTATS 2021. The sbibm suite.
#   https://github.com/sbi-benchmark/sbibm
#   neuralsbi ships four of its tasks as nsbi_task objects:
#   task_gaussian_linear(), task_two_moons(), task_slcp(), task_sir().
#   gaussian_linear has an analytic posterior, which lets us check that a
#   passing SBC really does correspond to a correct posterior. two_moons is
#   used below as the case where a fit is structurally unable to be right.
#
# Runtime: about 4 minutes on a laptop CPU.

library(neuralsbi)

has_torch <- requireNamespace("torch", quietly = TRUE) &&
  torch::torch_is_installed()
has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

# ---------------------------------------------------------------------------
# The task
# ---------------------------------------------------------------------------

task <- task_gaussian_linear(dim = 3L)
print(task)

# prior theta ~ N(0, 0.1 I); likelihood x | theta ~ N(theta, 0.1 I).
# reference_posterior(x_obs, n) returns exact draws.
prior <- task$prior
simulator <- task$simulator

# ---------------------------------------------------------------------------
# A well-trained fit
# ---------------------------------------------------------------------------

fit <- npe(prior, simulator, n_simulations = 4000,
           density_estimator = if (has_torch) "maf" else "linear_gaussian",
           seed = 1)

sbc_good <- sbc(fit, simulator, n_sbc = 300L, n_posterior_samples = 400L,
                seed = 1)
print(sbc_good)

# ---------------------------------------------------------------------------
# Reading the output
# ---------------------------------------------------------------------------
#
# $ranks is n_sbc x dim_theta: the rank of the true value among the posterior
# draws for that trial. Under calibration each column is uniform on
# 0..n_posterior_samples.

str(sbc_good$ranks)
print(head(sbc_good$ranks, 3))

# The p-values are a chi-square test of uniformity on binned ranks, computed by
# Monte Carlo so the asymptotic approximation never has to hold. Large is good.
# Treat them as a screen rather than a verdict: at a few hundred trials the
# test has modest power and is noisy, and one value under 0.05 out of three
# parameters is what you expect one time in twenty per parameter.
print(round(sbc_good$uniformity_pvalue, 3))

if (has_ggplot) {
  plot_sbc(sbc_good, param = 1L)
  # The dashed line is the expected count per bin; the dotted lines are a 99%
  # band from the binomial sampling noise of a finite number of trials.
}

# Coverage says the same thing in units people already use: of the trials, what
# fraction had the truth inside the nominal central interval? It is steadier
# than the p-value and it is what to report.
levels <- c(0.5, 0.8, 0.9, 0.95)
print(round(expected_coverage(sbc_good, levels), 3))
if (has_ggplot) plot_coverage(sbc_good)

# ---------------------------------------------------------------------------
# Does passing SBC mean the posterior is right?
# ---------------------------------------------------------------------------
#
# This task has an analytic posterior, so we can check both. Here they agree,
# which is the outcome to want.

set.seed(3)
theta_true <- sample_prior(prior, 1)
x_obs <- simulator(theta_true[1, ])

draws <- sample(posterior(fit, x_obs = x_obs), 4000)
ref <- task$reference_posterior(x_obs, n = 4000)

print(c2st(draws, ref, seed = 1))
print(round(rbind(npe = colMeans(draws), exact = colMeans(ref)), 3))
print(round(rbind(npe = apply(draws, 2, sd), exact = apply(ref, 2, sd)), 3))

# ---------------------------------------------------------------------------
# What a failing SBC looks like
# ---------------------------------------------------------------------------
#
# Starving a fit of simulations is the obvious way to break it, and on this
# task it does not work: the posterior stays close to correct because the
# model really is linear-Gaussian. To see SBC fire you need a fit that CANNOT
# be right, and the cleanest example is an estimator whose shape is wrong for
# the posterior it has to represent.
#
# two_moons has a curved, bimodal posterior. "linear_gaussian" is a
# closed-form conditional Gaussian: one mode, elliptical contours. It has no
# way to cover both crescents, and SBC says so.

tm <- task_two_moons()
fit_wrong <- npe(tm$prior, tm$simulator, n_simulations = 4000,
                 density_estimator = "linear_gaussian", seed = 1)
sbc_wrong <- sbc(fit_wrong, tm$simulator, n_sbc = 300L,
                 n_posterior_samples = 400L, seed = 1)
print(sbc_wrong)
print(round(expected_coverage(sbc_wrong, levels), 3))

if (has_ggplot) {
  plot_sbc(sbc_wrong, param = 1L)
  plot_coverage(sbc_wrong)
}

# Expect empirical coverage well below nominal at every level: a "90%
# interval" that contains the truth about 79% of the time. In the rank
# histogram that is the U shape, ranks piling up at both edges because the
# truth keeps landing outside a posterior that is too narrow.
#
# Read the shape as a diagnosis:
#
#   flat        calibrated
#   U           posterior too narrow, the failure that matters
#   arch        posterior too wide, conservative
#   sloped      posterior biased in that parameter
#   spike at 0  posterior systematically above the truth
#
# A U shape on a model whose estimator CAN represent the posterior is usually a
# simulation-budget problem. Raise n_simulations before touching architecture.

# ---------------------------------------------------------------------------
# The blind spot: SBC cannot see a posterior that ignores the data
# ---------------------------------------------------------------------------
#
# The averaging that makes SBC work also limits it. If the "posterior" is just
# the prior, the rank of a prior draw among prior draws is uniform, so SBC
# passes with nothing to report. This is not a hypothetical: an estimator that
# fails to learn anything from x lands exactly here.
#
# Construct it on purpose. Shuffle x against theta so no information survives,
# and fit as usual.

sims <- simulate_for_sbi(simulator, prior, n = 4000, seed = 2)
shuffled <- sims$x[base::sample.int(nrow(sims$x)), , drop = FALSE]
fit_blind <- npe(prior, theta = sims$theta, x = shuffled,
                 density_estimator = "linear_gaussian", seed = 1)

sbc_blind <- sbc(fit_blind, simulator, n_sbc = 300L,
                 n_posterior_samples = 400L, seed = 1)
print(sbc_blind)
print(round(expected_coverage(sbc_blind, levels), 3))

# Calibrated, by every SBC measure. And useless:
d_blind <- sample(posterior(fit_blind, x_obs = x_obs), 4000)
cat("\nposterior sd :", sprintf("%.3f", apply(d_blind, 2, sd)), "\n")
cat("prior sd     :", sprintf("%.3f", apply(sample_prior(prior, 4000), 2, sd)),
    "\n")
cat("exact post sd:", sprintf("%.3f", apply(ref, 2, sd)), "\n")
print(c2st(d_blind, ref, seed = 1))

# So run SBC to rule out miscalibration, and something else to rule out
# uselessness: a c2st() against a reference where one exists (10_c2st.R), a
# posterior predictive check where none does (13), and always a comparison of
# the posterior's width against the prior's.

# ---------------------------------------------------------------------------
# How many trials, how many draws
# ---------------------------------------------------------------------------
#
#   n_sbc                how many independent trials. This sets the power of
#                        the test. 100 is a smoke check, 300 to 1000 is a real
#                        one. Each trial costs one simulation plus one
#                        posterior sampling.
#   n_posterior_samples  the resolution of each rank. 100 is enough to see a
#                        strong U; a few hundred to resolve a mild slope.
#
# For an NPE fit the sampling is a forward pass and the whole thing is cheap.
# For NLE and NRE every trial is a separate MCMC run, so start small.

cheap <- sbc(fit, simulator, n_sbc = 100L, n_posterior_samples = 100L, seed = 2)
print(cheap)

# ---------------------------------------------------------------------------
# SBC on an NLE fit
# ---------------------------------------------------------------------------
#
# Same call. The MCMC controls reach posterior() through ..., which is how you
# keep 40 separate chains affordable.

fit_nle <- nle(prior, simulator, n_simulations = 2000,
               density_estimator = "linear_gaussian", seed = 1)
t0 <- Sys.time()
sbc_nle <- sbc(fit_nle, simulator, n_sbc = 40L, n_posterior_samples = 200L,
               n_chains = 10, warmup = 100, thin = 1, seed = 1)
cat(sprintf("\n40 NLE trials took %.0f s\n",
            as.numeric(Sys.time() - t0, units = "secs")))
print(sbc_nle)

# ---------------------------------------------------------------------------
# Two things that will trip you up
# ---------------------------------------------------------------------------
#
# 1. Leave `prior` alone. sbc() defaults to fit$prior, and that default is the
#    only choice that answers "is this fit calibrated". Passing a different
#    prior changes the question to how the fit behaves on parameters it was
#    never calibrated against, which is a fair thing to ask and is not SBC.
#
# 2. A trial whose posterior comes back short is an error, not a warning.
#    Ranks are binned against n_posterior_samples, so a trial scored on fewer
#    draws would read as miscalibration that is not there. It happens when a
#    bounded prior and a leaky estimator defeat rejection sampling; the fix is
#    a better fit, not a smaller n_posterior_samples.
