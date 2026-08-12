# 13_prior_posterior_predictive_checks.R -------------------------------------
#
# Two checks that bracket the fit. The prior predictive check runs BEFORE any
# inference: simulate from the prior and ask whether the model can produce data
# like yours at all. If it cannot, no estimator will save you, and you have
# learned that for the price of a few hundred simulations. The posterior
# predictive check runs after: simulate from the posterior and ask whether the
# fitted model reproduces the data it was conditioned on.
#
# Model source
#   Filipovic-Pierucci, A., Zarca, K. and Durand-Zaleski, I. heemod, vignette
#   "Time-homogeneous Markov model".
#   https://cran.r-project.org/web/packages/heemod/vignettes/c_homogeneous.html
#
#   The HIV therapy Markov model with four states: A (200 < CD4 < 500),
#   B (CD4 < 200), C (AIDS), D (death). The vignette's monotherapy transition
#   matrix, from Chancellor et al. (1997) via Briggs, Claxton and Sculpher,
#   "Decision Modelling for Health Economic Evaluation", is
#
#     A -> A .721   A -> B .202   A -> C .067   A -> D .010
#                   B -> B .581   B -> C .407   B -> D .012
#                                 C -> C .750   C -> D .250
#
#   The vignette takes those probabilities as known. Here we treat the five
#   free ones out of A and B as unknown and recover them from cross-sectional
#   surveys of the cohort, which is the situation you are in when the published
#   model has to be re-estimated against a different population. p_CD is fixed
#   at the vignette's 0.250.
#
# Runtime: about 3 minutes on a laptop CPU.

library(neuralsbi)

has_torch <- requireNamespace("torch", quietly = TRUE) &&
  torch::torch_is_installed()
has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

# ---------------------------------------------------------------------------
# The model
# ---------------------------------------------------------------------------

n_cycles <- 20            # 20 years
survey_at <- c(5, 10, 15, 20)
n_survey <- 200           # patients sampled at each survey
p_CD <- 0.250             # fixed at the vignette's value

markov_trace <- function(p_AB, p_AC, p_AD, p_BC, p_BD) {
  P <- matrix(0, 4, 4)
  P[1, ] <- c(1 - p_AB - p_AC - p_AD, p_AB, p_AC, p_AD)
  P[2, ] <- c(0, 1 - p_BC - p_BD, p_BC, p_BD)
  P[3, ] <- c(0, 0, 1 - p_CD, p_CD)
  P[4, ] <- c(0, 0, 0, 1)
  M <- matrix(0, n_cycles + 1L, 4L)
  M[1, ] <- c(1, 0, 0, 0)
  for (t in seq_len(n_cycles)) M[t + 1L, ] <- M[t, ] %*% P
  M
}

# One simulator call: run the cohort, then draw a survey of n_survey patients
# at each of four years and report the proportion in A, B and C. The proportion
# in D is one minus the rest, so reporting it would be redundant.
simulator <- function(p_AB, p_AC, p_AD, p_BC, p_BD) {
  M <- markov_trace(p_AB, p_AC, p_AD, p_BC, p_BD)[survey_at + 1L, , drop = FALSE]
  out <- numeric(0)
  for (i in seq_along(survey_at)) {
    counts <- rmultinom(1L, n_survey, pmax(M[i, ], 0))[, 1L]
    out <- c(out, counts[1:3] / n_survey)
  }
  stats::setNames(out, paste0(rep(c("A", "B", "C"), times = length(survey_at)),
                              "_y", rep(survey_at, each = 3)))
}

theta_true <- c(p_AB = 0.202, p_AC = 0.067, p_AD = 0.010,
                p_BC = 0.407, p_BD = 0.012)

set.seed(12)
x_obs <- do.call(simulator, as.list(theta_true))
print(round(x_obs, 3))

# ---------------------------------------------------------------------------
# 1. The prior predictive check
# ---------------------------------------------------------------------------
#
# Do this before fitting anything. simulate_for_sbi() draws from the prior and
# runs the simulator; the spread of the results is the prior predictive
# distribution. If the observation sits outside it, the model cannot produce
# your data under any parameter the prior allows, and the fit that follows will
# be a confident answer to the wrong question.
#
# Start with a prior that is too tight on purpose, to see the check fire.

prior_tight <- prior_uniform(
  low  = c(p_AB = 0.30, p_AC = 0.10, p_AD = 0.020, p_BC = 0.50, p_BD = 0.020),
  high = c(p_AB = 0.45, p_AC = 0.18, p_AD = 0.045, p_BC = 0.60, p_BD = 0.045)
)

prior_pred_check <- function(prior, label, n = 600) {
  sims <- simulate_for_sbi(simulator, prior, n = n, seed = 3)
  lo <- apply(sims$x, 2, quantile, 0.025)
  hi <- apply(sims$x, 2, quantile, 0.975)
  inside <- x_obs >= lo & x_obs <= hi
  cat(sprintf("\n%s: %d/%d outcomes inside the 95%% prior predictive band\n",
              label, sum(inside), length(inside)))
  print(data.frame(outcome = names(x_obs),
                   observed = round(x_obs, 3),
                   pp_lo = round(lo, 3),
                   pp_hi = round(hi, 3),
                   inside = inside, row.names = NULL))
  invisible(sims)
}

prior_pred_check(prior_tight, "tight prior")

# Several outcomes should fall outside. That is the check working: the prior
# rules out the region the data came from.

prior <- prior_uniform(
  low  = c(p_AB = 0.05, p_AC = 0.01, p_AD = 0.001, p_BC = 0.20, p_BD = 0.001),
  high = c(p_AB = 0.40, p_AC = 0.20, p_AD = 0.050, p_BC = 0.60, p_BD = 0.050)
)
sims <- prior_pred_check(prior, "wide prior")

# ---------------------------------------------------------------------------
# The plot version
# ---------------------------------------------------------------------------
#
# plot_posterior_predictive() is not fussy about where the draws came from.
# Pass prior predictive draws and it draws the prior predictive check, with the
# observation marked. The invisible return is the fraction of draws below the
# observation in each dimension, which is the tail probability. Values near 0
# or 1 are the outcomes the model struggles with.

if (has_ggplot) {
  q_prior <- plot_posterior_predictive(sims$x, x_obs)
  cat("\nprior predictive tail probabilities:\n")
  print(round(q_prior, 2))
}

# ---------------------------------------------------------------------------
# 2. Fit
# ---------------------------------------------------------------------------

fit <- npe(prior, theta = sims$theta, x = sims$x,
           density_estimator = if (has_torch) "maf" else "linear_gaussian",
           seed = 1)

# 600 simulations is a prior predictive budget, not an inference budget. Refit
# on a real one.
fit <- npe(prior, simulator, n_simulations = 4000,
           density_estimator = if (has_torch) "maf" else "linear_gaussian",
           seed = 1)

post <- posterior(fit, x_obs = x_obs)
draws <- sample(post, 4000)
print(summary(draws))
cat("\ntruth:", paste(names(theta_true), theta_true, sep = "=", collapse = "  "),
    "\n")

# p_AD and p_BD are the hard ones. Ten deaths per thousand per year barely move
# the state occupancies these surveys measure, so the posterior for them stays
# close to the prior. That is the design of the study talking, not the
# estimator failing, and the posterior predictive check below will not flag it,
# because a model with the wrong p_AD still reproduces this data.

# ---------------------------------------------------------------------------
# 3. The posterior predictive check
# ---------------------------------------------------------------------------

pred <- posterior_predictive(post, simulator, n = 800)
str(pred)

band <- apply(pred, 2, quantile, probs = c(0.025, 0.975))
inside <- x_obs >= band[1, ] & x_obs <= band[2, ]
cat(sprintf("\n%d/%d outcomes inside the 95%% posterior predictive band\n",
            sum(inside), length(inside)))

print(data.frame(outcome = names(x_obs),
                 observed = round(x_obs, 3),
                 pred_median = round(apply(pred, 2, median), 3),
                 pp_lo = round(band[1, ], 3),
                 pp_hi = round(band[2, ], 3),
                 row.names = NULL))

if (has_ggplot) {
  q_post <- plot_posterior_predictive(pred, x_obs)
  cat("\nposterior predictive tail probabilities:\n")
  print(round(q_post, 2))
}

# ---------------------------------------------------------------------------
# Prior against posterior, in predictive space
# ---------------------------------------------------------------------------
#
# The posterior predictive should be much tighter than the prior predictive in
# the outcomes the data speak to, and barely tighter in the ones they do not.
# The ratio of widths says which is which, and it is a more useful summary than
# any single parameter's posterior sd.

width <- function(m) apply(m, 2, function(z) diff(quantile(z, c(0.025, 0.975))))
print(data.frame(outcome = names(x_obs),
                 prior_width = round(width(sims$x), 3),
                 posterior_width = round(width(pred), 3),
                 ratio = round(width(pred) / width(sims$x), 2),
                 row.names = NULL))

# ---------------------------------------------------------------------------
# Conditioning on data the model cannot produce
# ---------------------------------------------------------------------------
#
# The failure mode worth rehearsing. Here is an observation from a cohort that
# progresses far faster than anything in the prior. The posterior still returns
# draws, and they still look like a posterior. The predictive check is what
# tells you not to believe them.

set.seed(13)
x_weird <- do.call(simulator, as.list(c(p_AB = 0.7, p_AC = 0.2, p_AD = 0.05,
                                        p_BC = 0.6, p_BD = 0.05)))
draws_weird <- sample(post, 2000, obs = x_weird)
print(summary(draws_weird))

pred_weird <- posterior_predictive(post, simulator, n = 500, x = x_weird)
band_w <- apply(pred_weird, 2, quantile, probs = c(0.025, 0.975))
cat(sprintf("\noutcomes inside the 95%% band: %d/%d\n",
            sum(x_weird >= band_w[1, ] & x_weird <= band_w[2, ]),
            length(x_weird)))

if (has_ggplot) plot_posterior_predictive(pred_weird, x_weird)

# Almost nothing should be inside. A posterior that looks confident and a
# predictive check that misses the data is the signature of an observation
# outside the training distribution, and the fix is the prior or the model, not
# the estimator.
