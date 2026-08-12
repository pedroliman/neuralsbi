# 05_nre.R -------------------------------------------------------------------
#
# Neural ratio estimation. npe() learns the posterior, nle() learns the
# likelihood, nre() learns neither: it trains a binary classifier to tell
# (theta, x) pairs drawn from the joint apart from mismatched pairs, and the
# optimal classifier's logit is log p(x | theta) / p(x). The posterior follows
# from Bayes' rule and is sampled with MCMC.
#
# Why bother: a density estimator has to describe the whole shape of a
# distribution, including places that make no difference to the posterior. A
# classifier only has to say which parameter explains the data better. That is
# the easier job when the data are counts, mixed types, or anything a flow
# handles badly.
#
# Model source
#   MRC Centre for Global Infectious Disease Analysis, mcstate package,
#   "Example: SIR model" vignette.
#   https://mrc-ide.github.io/mcstate/articles/sir_models.html
#
#   A discrete-time stochastic SIR with binomial transitions,
#
#     n_SI ~ Binomial(S, 1 - exp(-beta * I / N * dt))
#     n_IR ~ Binomial(I, 1 - exp(-gamma * dt))
#
#   with dt = 0.25 (four steps per observed day), N = 1010, I0 = 10, and the
#   vignette's parameter values beta = 0.2, gamma = 0.1. Daily incidence is
#   accumulated within the update and observed with Poisson noise, which is the
#   vignette's compare function. mcstate fits this with a particle filter and
#   pMCMC; here we fit it with a classifier instead.
#
# Runtime: about 3 minutes on a laptop CPU.

library(neuralsbi)

has_torch <- requireNamespace("torch", quietly = TRUE) &&
  torch::torch_is_installed()

# ---------------------------------------------------------------------------
# The simulator
# ---------------------------------------------------------------------------

N <- 1010
I0 <- 10
dt <- 0.25
n_days <- 40
steps_per_day <- round(1 / dt)

simulator <- function(beta, gamma) {
  S <- N - I0
  I <- I0
  cases <- numeric(n_days)
  for (d in seq_len(n_days)) {
    inc <- 0
    for (k in seq_len(steps_per_day)) {
      n_SI <- rbinom(1L, S, 1 - exp(-beta * I / N * dt))
      n_IR <- rbinom(1L, I, 1 - exp(-gamma * dt))
      S <- S - n_SI
      I <- I + n_SI - n_IR
      inc <- inc + n_SI
    }
    cases[d] <- rpois(1L, inc + 1e-6)
  }
  stats::setNames(cases, paste0("day", seq_len(n_days)))
}

prior <- prior_uniform(low  = c(beta = 0.05, gamma = 0.02),
                       high = c(beta = 0.60, gamma = 0.40))

# The observation: one epidemic at the vignette's parameter values.
set.seed(2)
theta_true <- c(beta = 0.2, gamma = 0.1)
x_obs <- simulator(theta_true[["beta"]], theta_true[["gamma"]])
print(x_obs)

# ---------------------------------------------------------------------------
# Fit the ratio estimator
# ---------------------------------------------------------------------------
#
# Four classifiers:
#   resnet    residual MLP. The default, and sbi's.
#   mlp       plain MLP.
#   linear    one linear layer on the raw inputs.
#   logistic  closed-form logistic regression on quadratic features. No torch,
#             and fast enough to be a useful baseline.
#
# num_atoms is the training objective's contrast count: per simulation the
# classifier scores the true parameter against num_atoms - 1 parameters
# borrowed from the rest of the minibatch, through a softmax. More atoms mean a
# harder problem and a sharper ratio, at a linear cost in forward passes.
# num_atoms = 10 is sbi's default and this one.

classifier <- if (has_torch) "resnet" else "logistic"

fit <- nre(prior, simulator, n_simulations = 4000,
           classifier = classifier, num_atoms = 10L,
           seed = 2024, verbose = TRUE)
print(fit)
invisible(summary(fit))

# ---------------------------------------------------------------------------
# log_ratio(): the learned ratio
# ---------------------------------------------------------------------------
#
# The ratio is calibrated up to a constant that depends on x, which is all a
# posterior needs. Differences at a fixed observation are meaningful; the
# absolute level is not.

cand <- rbind(truth = theta_true,
              slower = c(0.12, 0.10),
              faster = c(0.35, 0.10),
              short  = c(0.20, 0.30))
lr <- log_ratio(fit, cand, x_obs)
print(round(lr - max(lr), 2))

# The true parameters should score highest. Nothing guarantees it on a single
# noisy epidemic, but a fit where they score badly is a fit to distrust.

# Rows of x are independent observations, and the log ratio sums over them, the
# same way an NLE log-likelihood does. Two independent epidemics from the same
# parameters:
set.seed(3)
x_two <- rbind(x_obs, simulator(theta_true[["beta"]], theta_true[["gamma"]]))
print(round(log_ratio(fit, cand, x_two, sum_iid = FALSE), 1))

# ---------------------------------------------------------------------------
# Posterior
# ---------------------------------------------------------------------------
#
# As with nle(), the posterior needs MCMC. There is only one sampler here, the
# vectorized slice sampler, and no sampler argument: the "stan" route works by
# transpiling a density into Stan code, and a classifier is not a density.

t0 <- Sys.time()
post <- posterior(fit, x_obs, n_chains = 20, warmup = 200, thin = 2, seed = 7)
draws <- sample(post, 3000)
cat(sprintf("MCMC: %.1f s\n", as.numeric(Sys.time() - t0, units = "secs")))

print(post)
print(summary(draws))
cat("\ntruth: beta 0.2, gamma 0.1\n")

# R0 = beta / gamma, which the vignette's values put at 2.
R0 <- draws[, "beta"] / draws[, "gamma"]
cat(sprintf("R0: %.2f [%.2f, %.2f]\n", mean(R0),
            quantile(R0, 0.025), quantile(R0, 0.975)))

# ---------------------------------------------------------------------------
# Two epidemics are worth more than one
# ---------------------------------------------------------------------------

post2 <- posterior(fit, x_two, n_chains = 20, warmup = 200, thin = 2, seed = 7)
d2 <- sample(post2, 3000)
cat("\nposterior sd of beta, 1 epidemic :", sprintf("%.4f", sd(draws[, "beta"])),
    "\nposterior sd of beta, 2 epidemics:", sprintf("%.4f", sd(d2[, "beta"])),
    "\n")

# ---------------------------------------------------------------------------
# log_prob() and MAP
# ---------------------------------------------------------------------------
#
# Unnormalized, and it says so: the softmax objective is invariant to adding
# any function of x to the logit, so the ratio carries an unknown x-dependent
# constant.

print(round(log_prob(post, cand), 2))
print(round(map_estimate(post), 3))

# ---------------------------------------------------------------------------
# Saving
# ---------------------------------------------------------------------------

path <- tempfile(fileext = ".rds")
save_nre(fit, path)
fit2 <- load_nre(path)
cat("\nreloaded, same ratios:",
    isTRUE(all.equal(log_ratio(fit, cand, x_obs),
                     log_ratio(fit2, cand, x_obs))), "\n")
unlink(path)

# ---------------------------------------------------------------------------
# A cheap baseline worth running first
# ---------------------------------------------------------------------------
#
# "logistic" is a closed-form logistic regression on quadratic features. It
# needs no torch and fits in a second. If it already recovers the parameters,
# the problem is easy and the neural classifier is buying you little; if it
# fails badly where the resnet succeeds, that gap is the value of the network.

fit_log <- nre(prior, simulator, n_simulations = 4000,
               classifier = "logistic", seed = 2024)
post_log <- posterior(fit_log, x_obs, n_chains = 20, warmup = 200, seed = 7)
print(summary(sample(post_log, 2000)))
