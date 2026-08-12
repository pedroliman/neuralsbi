# 10_c2st.R ------------------------------------------------------------------
#
# The classifier two-sample test. Train a classifier to tell two sets of draws
# apart; report its cross-validated accuracy. Near 0.5, the two sets are
# indistinguishable. Near 1.0, they are not. This is the standard number the
# SBI literature reports when a reference posterior exists, and it is the
# metric the sbibm benchmark ranks methods by.
#
# Task source
#   Lueckmann, J.-M., Boelts, J., Greenberg, D., Goncalves, P. and Macke, J.
#   "Benchmarking Simulation-Based Inference", AISTATS 2021.
#   https://github.com/sbi-benchmark/sbibm
#   task_gaussian_linear() is sbibm's gaussian_linear: prior N(0, 0.1 I),
#   likelihood N(theta, 0.1 I). Conjugate, so reference_posterior() is exact
#   rather than a long MCMC run, which makes it the right task for calibrating
#   your intuition about what a given C2ST value means.
#
# Runtime: about 3 minutes on a laptop CPU.

library(neuralsbi)

has_torch <- requireNamespace("torch", quietly = TRUE) &&
  torch::torch_is_installed()

task <- task_gaussian_linear(dim = 4L)
prior <- task$prior
simulator <- task$simulator

set.seed(5)
theta_true <- sample_prior(prior, 1)
x_obs <- simulator(theta_true[1, ])

# ---------------------------------------------------------------------------
# Calibrating your reading of the number
# ---------------------------------------------------------------------------
#
# Before comparing anything real, see what c2st() returns on cases where you
# already know the answer.

ref_a <- task$reference_posterior(x_obs, n = 4000)
ref_b <- task$reference_posterior(x_obs, n = 4000)

# Two independent draws from the SAME distribution. This is the floor.
cat("same distribution, twice:\n")
print(c2st(ref_a, ref_b, seed = 1))

# A shift of a quarter of a standard deviation in one coordinate.
shifted <- ref_b
shifted[, 1] <- shifted[, 1] + 0.25 * sd(ref_b[, 1])
cat("\n0.25 sd shift in one of four coordinates:\n")
print(c2st(ref_a, shifted, seed = 1))

# A 50% wider posterior, same mean. The classifier here is linear, so it is
# close to blind to this: no hyperplane separates two clouds that share a
# centre. Python sbi uses an MLP and would see it. A 0.5 from this function is
# therefore the weaker of the two claims, and should be read next to the
# moments rather than instead of them.
wider <- sweep(sweep(ref_b, 2, colMeans(ref_b)), 2, 1.5, `*`)
wider <- sweep(wider, 2, colMeans(ref_b), `+`)
cat("\nsame mean, 50% wider:\n")
print(c2st(ref_a, wider, seed = 1))
print(round(rbind(ref = apply(ref_a, 2, sd), wider = apply(wider, 2, sd)), 3))

# The prior, against the posterior. This is the ceiling: two obviously
# different distributions.
cat("\nprior vs posterior:\n")
print(c2st(sample_prior(prior, 4000), ref_a, seed = 1))

# ---------------------------------------------------------------------------
# Scoring a real fit
# ---------------------------------------------------------------------------

fit <- npe(prior, simulator, n_simulations = 4000,
           density_estimator = if (has_torch) "maf" else "linear_gaussian",
           seed = 1)
draws <- sample(posterior(fit, x_obs = x_obs), 4000)

res <- c2st(draws, ref_a, seed = 1)
print(res)

# $fold_accuracy is the per-fold spread. A mean of 0.55 with folds ranging from
# 0.50 to 0.62 is a different claim from a mean of 0.55 with every fold at
# 0.55, and only the second one is a stable number.
print(round(res$fold_accuracy, 3))

# ---------------------------------------------------------------------------
# C2ST as a function of the simulation budget
# ---------------------------------------------------------------------------
#
# This is the curve the benchmark papers plot, and it is the most useful thing
# C2ST does in practice: it tells you whether more simulations would help
# before you spend them.

budgets <- c(250, 1000, 4000)
curve <- vapply(budgets, function(n) {
  f <- npe(prior, simulator, n_simulations = n,
           density_estimator = if (has_torch) "maf" else "linear_gaussian",
           seed = 1)
  d <- sample(posterior(f, x_obs = x_obs), 4000)
  c2st(d, ref_a, seed = 1)$accuracy
}, numeric(1))
print(data.frame(n_simulations = budgets, c2st = round(curve, 3)))

# Falling toward 0.5 means the budget is the binding constraint. Flat and well
# above 0.5 means it is not, and the estimator or the summary statistics are.

# ---------------------------------------------------------------------------
# Comparing estimators at a fixed budget
# ---------------------------------------------------------------------------
#
# Train all of them on the same simulations so the comparison is about the
# estimator and nothing else.

sims <- simulate_for_sbi(simulator, prior, n = 4000, seed = 7)
estimators <- if (has_torch) c("linear_gaussian", "mdn", "maf", "nsf") else
  "linear_gaussian"

scores <- vapply(estimators, function(de) {
  f <- npe(prior, theta = sims$theta, x = sims$x,
           density_estimator = de, seed = 1)
  d <- sample(posterior(f, x_obs = x_obs), 4000)
  c2st(d, ref_a, seed = 1)$accuracy
}, numeric(1))
print(round(sort(scores), 3))

# On a linear-Gaussian task the "linear_gaussian" estimator is exact, so it
# should sit at the floor and no flow should beat it. When one does, that is
# noise, not a finding.

# ---------------------------------------------------------------------------
# Arguments and the mistakes they catch
# ---------------------------------------------------------------------------
#
#   x, y     matrices of draws, rows = draws, columns = dimensions. Both must
#            be the same width. Row counts need not match: the larger set is
#            subsampled down to the smaller, because accuracy against
#            unbalanced classes is not a two-sample test. 8000 draws against
#            2000 identical ones scores 0.8 for a classifier that has learned
#            nothing except to answer with the bigger class.
#   n_folds  cross-validation folds. At least 2, and fewer than the number of
#            draws in the smaller set.
#   seed     the fold split and the subsampling are random.

# Unequal sizes are handled, not rejected.
print(c2st(draws[1:4000, ], ref_a[1:800, ], seed = 1)$accuracy)

# Different widths are an error, because two sets of draws of different
# quantities are not a two-sample test of anything.
res <- try(c2st(draws, ref_a[, 1:2], seed = 1), silent = TRUE)
cat(conditionMessage(attr(res, "condition")), "\n")

# ---------------------------------------------------------------------------
# When there is no reference posterior
# ---------------------------------------------------------------------------
#
# Most real models have none. Two things still work.
#
# 1. Agreement between two fits that should agree. If two independently seeded
#    estimators trained on the same budget give distinguishable posteriors, at
#    least one is wrong.

fit_b <- npe(prior, simulator, n_simulations = 4000,
             density_estimator = if (has_torch) "maf" else "linear_gaussian",
             seed = 99)
draws_b <- sample(posterior(fit_b, x_obs = x_obs), 4000)
cat("\ntwo seeds, same budget:\n")
print(c2st(draws, draws_b, seed = 1))

# 2. Long-run MCMC on a surrogate likelihood as the reference, for a task where
#    the true likelihood is out of reach but nle() gives you something to
#    sample properly. That is a comparison of samplers rather than of
#    posteriors, so treat the answer accordingly.
