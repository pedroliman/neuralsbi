# 01_basic_npe_example.R -----------------------------------------------------
#
# The shortest end-to-end neural posterior estimation run on real data: prior,
# simulator, fit, posterior, check.
#
# Model and data source
#   Grinsztajn, L., Semenova, E., Margossian, C. C. and Riou, J. "Bayesian
#   workflow for disease transmission modeling in Stan". Stan case study.
#   https://mc-stan.org/learn-stan/case-studies/boarding_school_case_study.html
#   Stan source: https://github.com/stan-dev/example-models/blob/master/knitr/disease_transmission/sir_negbin.stan
#
#   The case study fits an SIR model to an influenza A (H1N1) outbreak in a
#   British boarding school, January 1978. N = 763 boys, 512 of whom fell ill.
#   The data are the number confined to bed on each of 14 days, from the
#   outbreaks package (outbreaks::influenza_england_1978_school$in_bed),
#   originally Anonymous (1978), "Influenza in a boarding school",
#   British Medical Journal 1:587.
#
#   We keep the case study's model exactly: the same ODE, the same
#   negative-binomial observation model, the same priors, the same y0. The only
#   difference is the inference engine. Stan runs NUTS against the likelihood;
#   here the likelihood is never written down. We simulate from the model and
#   train a conditional density estimator on the simulations.
#
# Runtime: about 2 minutes on a laptop CPU.

library(neuralsbi)

# ---------------------------------------------------------------------------
# The data
# ---------------------------------------------------------------------------

cases <- c(3, 8, 26, 76, 225, 298, 258, 233, 189, 128, 68, 29, 14, 4)
n_days <- length(cases)
N <- 763          # boys at risk
i0 <- 1           # index case
y0 <- c(S = N - i0, I = i0, R = 0)

# ---------------------------------------------------------------------------
# The simulator
# ---------------------------------------------------------------------------
#
# dS/dt = -beta * I * S / N
# dI/dt =  beta * I * S / N - gamma * I
# dR/dt =  gamma * I
#
# observed cases[t] ~ NegBinomial2(mean = I(t), phi = 1 / phi_inv)
#
# The case study integrates with Stan's rk45. A fixed-step RK4 at dt = 0.25 is
# indistinguishable here and keeps the script dependency-free. One simulator
# call returns one 14-day outbreak.

sir_infected <- function(beta, gamma, dt = 0.25) {
  deriv <- function(y) {
    inf <- beta * y[2L] * y[1L] / N
    rec <- gamma * y[2L]
    c(-inf, inf - rec, rec)
  }
  y <- y0
  n_steps <- round(n_days / dt)
  out <- numeric(n_days)
  for (s in seq_len(n_steps)) {
    k1 <- deriv(y)
    k2 <- deriv(y + dt / 2 * k1)
    k3 <- deriv(y + dt / 2 * k2)
    k4 <- deriv(y + dt * k3)
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
  phi <- 1 / phi_inv
  # rnbinom's `mu` parameterization is Stan's neg_binomial_2.
  stats::setNames(rnbinom(n_days, size = phi, mu = mu + 1e-8),
                  paste0("day", seq_len(n_days)))
}

# ---------------------------------------------------------------------------
# The prior
# ---------------------------------------------------------------------------
#
# The case study's Stan program declares beta, gamma and phi_inv with
# <lower=0> and gives them
#
#   beta    ~ normal(2, 1)
#   gamma   ~ normal(0.4, 0.5)
#   phi_inv ~ exponential(5)
#
# so the first two are normals truncated at zero. prior_truncated() is the R
# spelling of Stan's T[0,]: it cuts the density at the bound and renormalizes by
# the mass it keeps. prior_independent() then multiplies the three marginals
# together, and the argument names become the parameter names.

prior <- prior_independent(
  beta    = prior_truncated(prior_normal(mean = 2, sd = 1), lower = 0),
  gamma   = prior_truncated(prior_normal(mean = 0.4, sd = 0.5), lower = 0),
  phi_inv = prior_exponential(rate = 5)
)
print(prior)

# ---------------------------------------------------------------------------
# Fit
# ---------------------------------------------------------------------------
#
# Masked autoregressive flow, the default, and the same estimator Python sbi
# defaults to. Without torch, swap in density_estimator = "linear_gaussian":
# it will run, but this posterior is not Gaussian in x and it will show.

estimator <- if (requireNamespace("torch", quietly = TRUE) &&
                 torch::torch_is_installed()) "maf" else "linear_gaussian"

fit <- npe(prior, simulator,
           n_simulations = 4000,
           density_estimator = estimator,
           seed = 2024, verbose = TRUE)
print(fit)

# ---------------------------------------------------------------------------
# Posterior
# ---------------------------------------------------------------------------

post <- posterior(fit, x_obs = cases)
draws <- sample(post, 5000)

print(summary(draws))

# The quantities the case study reports.
R0 <- draws[, "beta"] / draws[, "gamma"]
recovery_time <- 1 / draws[, "gamma"]

cat("\nR0            : ", sprintf("%.2f [%.2f, %.2f]",
    mean(R0), quantile(R0, 0.025), quantile(R0, 0.975)), "\n")
cat("recovery time : ", sprintf("%.2f [%.2f, %.2f] days",
    mean(recovery_time), quantile(recovery_time, 0.025),
    quantile(recovery_time, 0.975)), "\n")

# The case study's NUTS posterior, for reference (its Table of summaries):
# beta about 1.73, gamma about 0.54, R0 about 3.2, recovery time about 1.9 days.
# NPE trained on 4000 simulations should land in the same neighbourhood. It
# will not match to three decimals, and it is not supposed to: this is an
# approximate posterior from a finite simulation budget.

# ---------------------------------------------------------------------------
# Does the fit reproduce the outbreak?
# ---------------------------------------------------------------------------
#
# Push posterior draws back through the simulator and compare with the data.
# This is the check that catches a posterior which is confidently wrong.

pred <- posterior_predictive(post, simulator, n = 500)

band <- apply(pred, 2, quantile, probs = c(0.05, 0.5, 0.95))
print(data.frame(day = 1:n_days,
                 observed = cases,
                 pred_median = round(band[2, ]),
                 pred_lo = round(band[1, ]),
                 pred_hi = round(band[3, ])))

# Fraction of days where the observation falls inside the 90% predictive band.
inside <- mean(cases >= band[1, ] & cases <= band[3, ])
cat("\ndays inside the 90% predictive band:", sprintf("%.0f%%", 100 * inside),
    "\n")

# ---------------------------------------------------------------------------
# Plots, if ggplot2 and friends are installed
# ---------------------------------------------------------------------------

if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("GGally", quietly = TRUE) &&
    requireNamespace("ggdensity", quietly = TRUE)) {
  pairplot(draws)
}

# ---------------------------------------------------------------------------
# What was and was not amortized
# ---------------------------------------------------------------------------
#
# The fit is amortized over the whole prior: conditioning on a different
# 14-day outbreak costs a forward pass, not a refit. Here is the same estimator
# conditioned on a hypothetical milder outbreak.

milder <- round(cases * 0.4)
print(summary(sample(post, 2000, obs = milder)))
