# 02_priors.R ----------------------------------------------------------------
#
# Everything neuralsbi asks of a prior, and the ways to build one:
# prior_uniform(), prior_normal(), the named families (prior_lognormal(),
# prior_beta(), prior_gamma() and the rest), prior_independent(),
# prior_truncated(), and prior_custom() for anything left over. Plus
# sample_prior(), within_support(), and how parameter names travel from the
# prior all the way to the axis labels of a diagnostic plot.
#
# Model source
#   Alarid-Escudero, F., Krijkamp, E., Enns, E. A., Yang, A., Hunink, M. G. M.,
#   Pechlivanoglou, P. and Jalal, H. "An Introductory Tutorial on Cohort
#   State-Transition Models in R for Cost-Effectiveness Analysis".
#   Medical Decision Making 43(1), 2023. Code:
#   https://github.com/DARTH-git/cohort-modeling-tutorial-intro/blob/main/analysis/cSTM_time_indep.R
#
#   The parameters and their calibration ranges are the Sick-Sicker model's, as
#   calibrated in the DARTH group's darthpack:
#   https://github.com/DARTH-git/darthpack
#   That package calibrates p_S1S2, hr_S1 and hr_S2 against three
#   epidemiological targets, over the ranges used below.
#
# Runtime: a few seconds.

library(neuralsbi)

# ---------------------------------------------------------------------------
# What a prior has to provide
# ---------------------------------------------------------------------------
#
# An nsbi_prior knows how to draw samples and (for the methods that need it)
# how to evaluate its own log-density. Bounded priors also carry lower/upper
# limits, which the posterior uses to reject out-of-support draws.

# ---------------------------------------------------------------------------
# 1. prior_uniform(): the box prior
# ---------------------------------------------------------------------------
#
# The three calibrated Sick-Sicker parameters and darthpack's search ranges.

prior <- prior_uniform(
  low  = c(p_S1S2 = 0.01, hr_S1 = 1.0, hr_S2 =  5),
  high = c(p_S1S2 = 0.50, hr_S1 = 4.5, hr_S2 = 15)
)
print(prior)

# Naming `low` is worth doing. Those names become the column names of every
# parameter matrix, posterior sample and SBC rank matrix downstream, and they
# decide how the simulator gets called (see 03_npe.R).

theta <- sample_prior(prior, 5)
print(theta)

# The log-density is a plain function of a matrix, one value per row.
print(prior$log_prob(theta))

# A box prior has hard edges, and within_support() is what tests them.
edge <- rbind(c(0.25, 3, 10),      # inside
              c(0.25, 3, 20),      # hr_S2 too high
              c(0.00, 3, 10))      # p_S1S2 too low
print(within_support(prior, edge))

# ---------------------------------------------------------------------------
# 2. prior_normal(): independent normals
# ---------------------------------------------------------------------------
#
# `sd` recycles, so one number covers every parameter.

prior_n <- prior_normal(mean = c(log_p_S1S2 = log(0.105),
                                 log_hr_S1  = log(3),
                                 log_hr_S2  = log(10)),
                        sd = 0.4)
print(prior_n)
print(round(sample_prior(prior_n, 3), 3))

# A normal prior is unbounded, so within_support() is TRUE everywhere. Nothing
# will be rejected, and nothing needs to be.
print(within_support(prior_n, c(-99, 99, 0)))

# ---------------------------------------------------------------------------
# 3. Named families: one marginal per parameter
# ---------------------------------------------------------------------------
#
# The natural prior for these three parameters is not a box. p_S1S2 is a
# probability, and the two hazard ratios are positive and right-skewed. A beta
# and two log-normals say that directly. Each family constructor takes Stan's
# argument names, and each is vectorized, so one call covers both hazard
# ratios.

prior_hr <- prior_lognormal(meanlog = c(hr_S1 = log(3), hr_S2 = log(10)),
                            sdlog = c(0.3, 0.25))
print(prior_hr)
print(round(sample_prior(prior_hr, 3), 3))

# The support comes from the family, not from you. A log-normal is positive
# and a beta lives on [0, 1], and those bounds are what the posterior's
# leakage correction rejects and renormalizes against.
print(prior_hr$lower)
print(prior_beta(2, 15)$upper)

# ---------------------------------------------------------------------------
# 4. prior_independent(): a different family per parameter
# ---------------------------------------------------------------------------
#
# Stan writes a joint prior as one sampling statement per parameter and lets
# the product take care of itself. prior_independent() is that, as an object.
# Name the arguments and the names label every parameter matrix, posterior
# sample and diagnostic plot downstream.

prior_f <- prior_independent(
  p_S1S2 = prior_beta(2, 15),
  hr_S1  = prior_lognormal(log(3), 0.3),
  hr_S2  = prior_lognormal(log(10), 0.25)
)
print(prior_f)
print(round(sample_prior(prior_f, 5), 3))

# The log-density is the sum of the marginals, one value per row.
print(prior_f$log_prob(sample_prior(prior_f, 3)))

# ---------------------------------------------------------------------------
# 5. prior_truncated(): Stan's T[lower, upper]
# ---------------------------------------------------------------------------
#
# darthpack searches hr_S2 between 5 and 15. Truncating the log-normal to that
# interval renormalizes the density by the mass it keeps, so log_prob() stays a
# proper log density instead of the untruncated one shifted by an unknown
# constant. That matters here in a way it does not in Stan, where the constant
# drops out: the prior density is summed with a learned likelihood in nle() and
# compared against reference draws in the diagnostics.

prior_t <- prior_truncated(prior_lognormal(log(10), 0.25),
                           lower = 5, upper = 15)
print(prior_t)
print(range(sample_prior(prior_t, 1000)))
print(integrate(function(z) exp(prior_t$log_prob(z)), 5, 15)$value)

# A half-normal is the same idea with one bound, and comes ready-made.
print(prior_half_normal(sd = 2))

# ---------------------------------------------------------------------------
# 6. prior_custom(): anything else
# ---------------------------------------------------------------------------
#
# You supply the sampler and (optionally) the log-density. Both are checked at
# construction with one probe call, so a shape mistake is caught here rather
# than after the simulation budget is spent.
#
# Prefer a named family wherever one fits. A family carries its own support
# bounds, survives prior_truncated() with an exact normalizing constant, and is
# what stan_code() writes out as a sampling statement; prior_custom() is
# arbitrary R code and does none of those. What it buys you is dependence
# between parameters, which no product of marginals can express. Here the
# Sicker hazard ratio is forced above the Sick one.

prior_c <- prior_custom(
  sample_fn = function(n) {
    hr_S1 <- stats::rlnorm(n, log(3), 0.3)
    cbind(hr_S1 = hr_S1, hr_S2 = hr_S1 * (1 + stats::rexp(n)))
  },
  log_prob_fn = function(theta) {
    hr_S1 <- theta[, 1]
    hr_S2 <- theta[, 2]
    lp <- stats::dlnorm(hr_S1, log(3), 0.3, log = TRUE) +
      stats::dexp(hr_S2 / hr_S1 - 1, log = TRUE) - log(hr_S1)
    ifelse(hr_S2 > hr_S1, lp, -Inf)
  },
  dim = 2,
  lower = c(0, 0),
  param_names = c("hr_S1", "hr_S2")
)
print(prior_c)
draws_c <- sample_prior(prior_c, 5)
print(round(draws_c, 3))
print(all(draws_c[, "hr_S2"] > draws_c[, "hr_S1"]))

# The three mistakes prior_custom() catches immediately. Uncomment to see the
# messages.
#
# prior_custom(sample_fn = function(n) rnorm(n), dim = 2)
#   -> sample_fn(2) must return a 2 x 2 matrix
#
# prior_custom(sample_fn = function(n) cbind(rnorm(n), rnorm(n)),
#              log_prob_fn = function(theta) 0, dim = 2)
#   -> log_prob_fn must return one log-density per row
#
# prior_custom(sample_fn = function(n) cbind(rnorm(n), rnorm(n)),
#              dim = 2, lower = c(0, 0, 0))
#   -> `lower` must have one bound per parameter

# log_prob_fn is optional. Without it, sampling and NPE still work; MCMC-based
# posteriors (nle(), nre()) and anything that needs a prior density do not.
prior_nolp <- prior_custom(
  sample_fn = function(n) cbind(stats::rbeta(n, 2, 15)),
  dim = 1, lower = 0, upper = 1
)
print(prior_nolp$log_prob(matrix(0.1)))   # NA, by design

# ---------------------------------------------------------------------------
# 7. Names flow through everything
# ---------------------------------------------------------------------------
#
# Here is the Sick-Sicker cohort model from the DARTH tutorial, cut down to the
# three calibrated parameters, so we can watch the names travel.

n_cycles <- 75          # ages 25 to 100, annual cycles
r_HD  <- 0.002          # background mortality rate when Healthy
r_HS1 <- 0.15           # Healthy -> Sick
r_S1H <- 0.5            # Sick -> Healthy

rate_to_prob <- function(r, t = 1) 1 - exp(-r * t)

sick_sicker_trace <- function(p_S1S2, hr_S1, hr_S2) {
  p_HS1 <- rate_to_prob(r_HS1)
  p_S1H <- rate_to_prob(r_S1H)
  p_HD  <- rate_to_prob(r_HD)
  p_S1D <- rate_to_prob(r_HD * hr_S1)
  p_S2D <- rate_to_prob(r_HD * hr_S2)

  P <- matrix(0, 4, 4, dimnames = list(c("H", "S1", "S2", "D"),
                                       c("H", "S1", "S2", "D")))
  P["H",  "H"]  <- (1 - p_HD) * (1 - p_HS1)
  P["H",  "S1"] <- (1 - p_HD) * p_HS1
  P["H",  "D"]  <- p_HD
  P["S1", "H"]  <- (1 - p_S1D) * p_S1H
  P["S1", "S1"] <- (1 - p_S1D) * (1 - (p_S1H + p_S1S2))
  P["S1", "S2"] <- (1 - p_S1D) * p_S1S2
  P["S1", "D"]  <- p_S1D
  P["S2", "S2"] <- 1 - p_S2D
  P["S2", "D"]  <- p_S2D
  P["D",  "D"]  <- 1

  M <- matrix(0, n_cycles + 1L, 4L,
              dimnames = list(NULL, c("H", "S1", "S2", "D")))
  M[1, ] <- c(1, 0, 0, 0)
  for (t in seq_len(n_cycles)) M[t + 1L, ] <- M[t, ] %*% P
  M
}

# Survival at four ages, as one summary of the cohort.
target_cycles <- c(10, 25, 40, 55)

simulator <- function(p_S1S2, hr_S1, hr_S2) {
  M <- sick_sicker_trace(p_S1S2, hr_S1, hr_S2)
  stats::setNames(1 - M[target_cycles + 1L, "D"],
                  paste0("surv_age", 25 + target_cycles))
}

# simulate_for_sbi() draws from the prior and runs the simulator once per draw.
sims <- simulate_for_sbi(simulator, prior, n = 200, seed = 1)
str(sims)

# Both matrices are named: theta from the prior, x from the simulator's output.
print(head(round(sims$theta, 3), 3))
print(head(round(sims$x, 4), 3))

# ---------------------------------------------------------------------------
# 8. Bounded priors and posterior leakage
# ---------------------------------------------------------------------------
#
# A density estimator trained in unconstrained space puts some mass outside a
# bounded prior's support. The posterior handles that by rejection sampling and
# by renormalizing log_prob() by the estimated acceptance probability. That
# machinery is driven entirely by the prior's lower/upper. prior_uniform()
# carries them from its box, a named family from its own support, and
# prior_custom() only if you pass them; prior_normal() is unbounded and carries
# none.

fit <- npe(prior, simulator, theta = sims$theta, x = sims$x,
           density_estimator = "linear_gaussian", seed = 1)
post <- posterior(fit, x_obs = simulator(0.105, 3, 10))
draws <- sample(post, 1000)

# Every returned draw is inside the box, by construction.
cat("all draws in support:", all(within_support(prior, draws)), "\n")
# The acceptance rate that got them there is recorded on the draws.
print(draws)
