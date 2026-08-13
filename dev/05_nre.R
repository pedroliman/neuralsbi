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
# Runtime: about 4 minutes on a laptop CPU.

library(neuralsbi)

has_torch <- requireNamespace("torch", quietly = TRUE) &&
  torch::torch_is_installed()

# ---------------------------------------------------------------------------
# The simulator
# ---------------------------------------------------------------------------

N <- 1010
I0 <- 10
dt <- 0.25
n_days <- 100
steps_per_day <- round(1 / dt)

epidemic <- function(beta, gamma) {
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
  cases
}

# The simulator returns log(1 + cases), not the counts. Two reasons, one solid
# and one measured.
#
# The solid one: for a negative-binomial or Poisson count the noise scale grows
# with the mean, so a single per-column standard deviation is right at one
# point on the epidemic curve and wrong everywhere else. On the log scale the
# noise is roughly constant. And in the growth phase the parameter of interest
# is the slope of log(cases), which is an additive feature of the log series
# and a multiplicative one of the raw series. BayesFlow applies log1p to count
# data for the same reason and turns its own standardization off afterwards
# (examples/SIR_Posterior_Estimation.ipynb).
#
# The measured one, and it is smaller than the argument above would suggest.
# Fitting this model three times per representation, at seeds 2024, 7 and 99,
# against a truth of beta = 0.20:
#
#   raw counts  beta 0.101, 0.128, 0.129   posterior sd about 0.05
#   log1p       beta 0.141, 0.153, 0.212   posterior sd about 0.075
#
# The log fits are closer in all three, but mostly by being wider rather than
# better centred: the truth falls inside the 95% interval 3 times out of 3 on
# the log scale and 2 out of 3 on the raw one. Both are biased low. The raw
# fits are confidently wrong, which is the worse failure of the two, and that
# is the honest case for the transform here.
#
# Two things this is NOT. It is not that z-scoring destroys the signal: a
# k-nearest-neighbour lookup for beta in the standardized raw series does as
# well as in the log one (RMSE 0.065 against 0.065, prior sd 0.151). And it
# does not carry over to npe() on this model, where the transform makes no
# difference and the fit is simply short of simulations in 100 dimensions.
simulator <- function(beta, gamma) {
  stats::setNames(log1p(epidemic(beta, gamma)), paste0("day", seq_len(n_days)))
}

prior <- prior_uniform(low  = c(beta = 0.05, gamma = 0.02),
                       high = c(beta = 0.60, gamma = 0.40))

# The observation: one epidemic at the vignette's parameter values.
set.seed(2)
theta_true <- c(beta = 0.2, gamma = 0.1)
cases_obs <- epidemic(theta_true[["beta"]], theta_true[["gamma"]])
x_obs <- stats::setNames(log1p(cases_obs), paste0("day", seq_len(n_days)))
cat("total cases:", sum(cases_obs), " peak on day:", which.max(cases_obs), "\n")
print(cases_obs)

# ---------------------------------------------------------------------------
# Fit the ratio estimator
# ---------------------------------------------------------------------------
#
# Four classifiers:
#   resnet    residual MLP. The default, and sbi's.
#   mlp       plain MLP.
#   linear    one linear layer on the raw inputs.
#   logistic  closed-form logistic regression on quadratic features. No torch.
#             Suited to low-dimensional x; see the note at the end.
#
# num_atoms is the training objective's contrast count: per simulation the
# classifier scores the true parameter against num_atoms - 1 parameters
# borrowed from the rest of the minibatch, through a softmax. More atoms mean a
# harder problem and a sharper ratio, at a linear cost in forward passes.
# num_atoms = 10 is sbi's default and this one.

classifier <- if (has_torch) "resnet" else "logistic"

# Simulate once and reuse, so the classifiers compared below all see the same
# epidemics. max_epochs and patience are trimmed a little for runtime; the
# defaults are 2000 and 20.
sims <- simulate_for_sbi(simulator, prior, n = 4000, seed = 11)
cat("simulated", nrow(sims$x), "epidemics,", sims$n_dropped, "dropped\n")

fit <- nre(prior, theta = sims$theta, x = sims$x,
           classifier = classifier, num_atoms = 10L,
           max_epochs = 300L, patience = 15L,
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
# The other classifiers
# ---------------------------------------------------------------------------
#
# "mlp" is a plain multilayer perceptron and "linear" is a single linear layer
# on the raw (theta, x). The linear one is the useful control: if it already
# recovers the parameters, the discrimination problem is easy and the residual
# network is buying you little.

if (has_torch) {
  for (cl in c("mlp", "linear")) {
    f <- nre(prior, theta = sims$theta, x = sims$x, classifier = cl,
             max_epochs = 300L, patience = 15L, seed = 2024)
    d <- sample(posterior(f, x_obs, n_chains = 20, warmup = 150, seed = 7),
                2000)
    cat(sprintf("%-8s beta %.3f (sd %.3f)   gamma %.3f (sd %.3f)\n",
                cl, mean(d[, "beta"]), sd(d[, "beta"]),
                mean(d[, "gamma"]), sd(d[, "gamma"])))
  }
}

# The fourth option, "logistic", is a closed-form logistic regression on
# quadratic features and needs no torch. It is the fallback when torch is
# unavailable, and it is meant for low-dimensional x: the feature expansion is
# quadratic in dim_theta + dim_x, so on the 100-day series here it would build
# more than five thousand features per pair. Aggregate to weekly counts first,
# or use it on a model with a handful of summary statistics.
