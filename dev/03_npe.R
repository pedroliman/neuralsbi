# 03_npe.R -------------------------------------------------------------------
#
# Neural posterior estimation in depth: the four density estimators, the
# training controls that matter, saving and reloading a fit, and everything a
# posterior object can do (sample, log_prob, map_estimate, summary,
# as.data.frame). Ends with npe_sequential(), which trades amortization for a
# sharper posterior at one observation.
#
# Model source
#   Alarid-Escudero, F., Krijkamp, E., Enns, E. A., Yang, A., Hunink, M. G. M.,
#   Pechlivanoglou, P. and Jalal, H. "An Introductory Tutorial on Cohort
#   State-Transition Models in R for Cost-Effectiveness Analysis".
#   Medical Decision Making 43(1), 2023.
#   https://github.com/DARTH-git/cohort-modeling-tutorial-intro/blob/main/analysis/cSTM_time_indep.R
#
#   The Sick-Sicker model: four states (Healthy, Sick, Sicker, Dead), annual
#   cycles from age 25 to 100. Fixed inputs are the tutorial's (r_HD = 0.002,
#   r_HS1 = 0.15, r_S1H = 0.5). The three parameters calibrated here, and the
#   ranges they are searched over, are the ones the DARTH group's darthpack
#   calibrates: p_S1S2, hr_S1, hr_S2.
#   https://github.com/DARTH-git/darthpack
#
#   darthpack calibrates against three epidemiological targets: overall
#   survival, prevalence of disease, and the proportion of the sick who are in
#   the Sicker state. We use the same three, observed at four ages in a cohort
#   of 1000, with the binomial sampling noise a real survey would carry.
#
# Runtime: about 3 minutes on a laptop CPU.

library(neuralsbi)

has_torch <- requireNamespace("torch", quietly = TRUE) &&
  torch::torch_is_installed()

# ---------------------------------------------------------------------------
# The model
# ---------------------------------------------------------------------------

n_cycles <- 75
r_HD  <- 0.002
r_HS1 <- 0.15
r_S1H <- 0.5
n_cohort <- 1000                 # survey size behind each target
target_cycles <- c(10, 25, 40, 55)

rate_to_prob <- function(r, t = 1) 1 - exp(-r * t)

sick_sicker_trace <- function(p_S1S2, hr_S1, hr_S2) {
  p_HS1 <- rate_to_prob(r_HS1)
  p_S1H <- rate_to_prob(r_S1H)
  p_HD  <- rate_to_prob(r_HD)
  p_S1D <- rate_to_prob(r_HD * hr_S1)
  p_S2D <- rate_to_prob(r_HD * hr_S2)

  P <- matrix(0, 4, 4)
  P[1, ] <- c((1 - p_HD) * (1 - p_HS1), (1 - p_HD) * p_HS1, 0, p_HD)
  P[2, ] <- c((1 - p_S1D) * p_S1H,
              (1 - p_S1D) * (1 - (p_S1H + p_S1S2)),
              (1 - p_S1D) * p_S1S2, p_S1D)
  P[3, ] <- c(0, 0, 1 - p_S2D, p_S2D)
  P[4, ] <- c(0, 0, 0, 1)

  M <- matrix(0, n_cycles + 1L, 4L)
  M[1, ] <- c(1, 0, 0, 0)
  for (t in seq_len(n_cycles)) M[t + 1L, ] <- M[t, ] %*% P
  colnames(M) <- c("H", "S1", "S2", "D")
  M
}

# One simulator call: run the cohort, read the three targets at four ages, and
# add the binomial noise of a finite survey. Twelve numbers out.
simulator <- function(p_S1S2, hr_S1, hr_S2) {
  M <- sick_sicker_trace(p_S1S2, hr_S1, hr_S2)[target_cycles + 1L, , drop = FALSE]
  surv  <- 1 - M[, "D"]
  prev  <- (M[, "S1"] + M[, "S2"]) / surv
  props <- M[, "S2"] / (M[, "S1"] + M[, "S2"])

  noisy <- function(p, n) rbinom(length(p), n, pmin(pmax(p, 0), 1)) / n
  ages <- 25 + target_cycles
  c(stats::setNames(noisy(surv,  n_cohort), paste0("surv", ages)),
    stats::setNames(noisy(prev,  n_cohort), paste0("prev", ages)),
    stats::setNames(noisy(props, n_cohort), paste0("props", ages)))
}

prior <- prior_uniform(
  low  = c(p_S1S2 = 0.01, hr_S1 = 1.0, hr_S2 =  5),
  high = c(p_S1S2 = 0.50, hr_S1 = 4.5, hr_S2 = 15)
)

# The tutorial's own values, which we treat as the truth behind the data.
theta_true <- c(p_S1S2 = rate_to_prob(0.105), hr_S1 = 3, hr_S2 = 10)
set.seed(7)
x_obs <- do.call(simulator, as.list(theta_true))
print(round(x_obs, 3))

# ---------------------------------------------------------------------------
# How the simulator gets called
# ---------------------------------------------------------------------------
#
# The prior's parameter names are p_S1S2, hr_S1, hr_S2 and the simulator's
# formals are exactly those, so each parameter arrives as its own named
# argument. Any other simulator receives the whole parameter vector as its
# first argument. Both work; the named form is easier to read.

vec_simulator <- function(theta) simulator(theta[1], theta[2], theta[3])
str(simulate_for_sbi(vec_simulator, prior, n = 3, seed = 1)$x)

# ---------------------------------------------------------------------------
# Simulate once, fit many times
# ---------------------------------------------------------------------------
#
# npe() will call the simulator for you, but when comparing estimators you want
# them all trained on the same simulations. simulate_for_sbi() gives you the
# (theta, x) pair to pass in directly.

sims <- simulate_for_sbi(simulator, prior, n = 4000, seed = 42)
cat("simulated:", nrow(sims$x), "draws,", sims$n_dropped, "dropped\n")

# ---------------------------------------------------------------------------
# The four density estimators
# ---------------------------------------------------------------------------
#
#   linear_gaussian  closed-form conditional Gaussian. No torch. Exact when the
#                    model is linear-Gaussian, a useful sanity baseline when it
#                    is not.
#   mdn              mixture density network: an MLP emitting the weights,
#                    means and full covariances of a Gaussian mixture.
#   maf              masked autoregressive flow. The default, and sbi's.
#   nsf              autoregressive rational-quadratic spline flow. The most
#                    flexible and the slowest.

estimators <- if (has_torch) c("linear_gaussian", "mdn", "maf", "nsf") else
  "linear_gaussian"

fits <- lapply(estimators, function(de) {
  cat("\n---", de, "---\n")
  t0 <- Sys.time()
  f <- npe(prior, theta = sims$theta, x = sims$x,
           density_estimator = de, seed = 1)
  cat(sprintf("  %.1f s\n", as.numeric(Sys.time() - t0, units = "secs")))
  f
})
names(fits) <- estimators

# summary() prints the fit and returns its training metadata invisibly.
for (de in estimators) {
  cat("\n")
  info <- summary(fits[[de]])
}

# Posterior means, side by side against the truth.
means <- t(vapply(fits, function(f) {
  colMeans(sample(posterior(f, x_obs = x_obs), 2000))
}, numeric(3)))
print(round(rbind(means, truth = theta_true), 3))

# The flow estimators should recover p_S1S2 and hr_S2 well. hr_S1 is the hard
# one: with r_HD = 0.002 the Sick state's mortality barely moves survival, so
# the data carry little information about it and the posterior stays close to
# the prior. That is the model talking, not the estimator failing, and 09_sbc.R
# shows how to tell those two apart.

fit <- fits[[if (has_torch) "maf" else "linear_gaussian"]]

# ---------------------------------------------------------------------------
# Training controls
# ---------------------------------------------------------------------------
#
# The defaults match Python sbi: batch_size = 200, lr = 5e-4,
# validation_fraction = 0.1, patience = 20, and a max_epochs cap that early
# stopping normally reaches first. The two worth knowing:
#
#   n_restarts     train k independently initialized networks, keep the best
#                  validation loss. Guards against a bad initialization and
#                  against MDN mode collapse.
#   n_transforms   depth of a MAF or NSF. More transforms, more flexibility,
#                  more time.

if (has_torch) {
  quick <- npe(prior, theta = sims$theta, x = sims$x,
               density_estimator = "maf",
               n_transforms = 3L, hidden = c(30L, 30L),
               batch_size = 500L, patience = 10L,
               seed = 1)
  cat("\nsmaller, faster MAF:\n")
  invisible(summary(quick))
}

# ---------------------------------------------------------------------------
# What a posterior can do
# ---------------------------------------------------------------------------

post <- posterior(fit, x_obs = x_obs)
print(post)

# sample(): draws in the original parameter units, with names attached.
draws <- sample(post, 4000)
print(summary(draws))

# sample_posterior() is the same thing under a name that does not collide with
# base::sample, for use inside packages.
draws2 <- sample_posterior(post, n = 500)

# as.data.frame(): straight into a data-frame workflow.
df <- as.data.frame(draws)
str(df)

# log_prob(): the posterior log-density at parameter values you choose,
# renormalized for the mass the estimator leaks outside the prior box.
lp <- log_prob(post, rbind(theta_true, c(0.30, 3, 10), c(0.90, 3, 10)))
print(round(lp, 2))
# The third row is outside the prior support, so it is -Inf.

# map_estimate(): the highest-density parameter vector, found by optimizing
# log_prob from the best of n_init posterior draws.
print(round(map_estimate(post), 3))

# summary() on the posterior samples n draws and summarizes them.
print(summary(post, n = 2000))

# ---------------------------------------------------------------------------
# Amortization
# ---------------------------------------------------------------------------
#
# The fit is trained once over the whole prior. Conditioning on a different
# observation costs one forward pass, so a fit can serve many data sets. Here
# are three cohorts with different underlying progression rates.

for (p in c(0.05, 0.15, 0.30)) {
  set.seed(11)
  xi <- simulator(p, 3, 10)
  di <- sample(post, 2000, obs = xi)
  cat(sprintf("true p_S1S2 = %.2f -> posterior mean %.3f [%.3f, %.3f]\n",
              p, mean(di[, 1]),
              quantile(di[, 1], 0.025), quantile(di[, 1], 0.975)))
}

# ---------------------------------------------------------------------------
# Saving a fit
# ---------------------------------------------------------------------------
#
# A torch-backed fit holds an external pointer. saveRDS() writes the pointer,
# not the network: the file reloads and prints fine, then fails on first use.
# save_npe() writes the weights alongside everything else.

path <- tempfile(fileext = ".rds")
save_npe(fit, path)
fit2 <- load_npe(path)
print(round(colMeans(sample(posterior(fit2, x_obs = x_obs), 1000)), 3))
unlink(path)

# ---------------------------------------------------------------------------
# npe_sequential(): spend the budget where the posterior is
# ---------------------------------------------------------------------------
#
# Single-round npe() spreads its simulations over the whole prior. When only
# one observation matters, most of those land where the posterior never goes.
# npe_sequential() truncates the prior to the current posterior's
# high-probability region after each round and simulates there instead.
#
# The price: the result is no longer amortized. It is trustworthy at x_obs and
# nowhere else.

snpe <- npe_sequential(prior, simulator, x_obs = x_obs,
                       n_rounds = 2L, n_simulations = 1500L,
                       density_estimator = if (has_torch) "maf" else
                         "linear_gaussian",
                       seed = 3, verbose = TRUE)
print(snpe)

post_seq <- posterior(snpe, x_obs = x_obs)
print(summary(post_seq, n = 2000))

# Per-round simulation counts, acceptance rates and truncation thresholds.
print(do.call(rbind, lapply(seq_along(snpe$rounds), function(i) {
  r <- snpe$rounds[[i]]
  data.frame(round = i, n_new = r$n_new,
             acceptance = round(r$acceptance, 3),
             threshold = round(r$threshold, 2))
})))
