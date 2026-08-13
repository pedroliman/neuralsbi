# 15_npe_vs_pomp.R -----------------------------------------------------------
#
# The same partially observed Markov process, fit two ways: by iterated
# filtering in pomp, and by neural posterior estimation here. They are not
# competitors so much as different bargains, and this script is about what each
# one buys.
#
#   pomp   estimates the likelihood with a particle filter and maximizes it
#          with iterated filtering. Every likelihood evaluation costs a full
#          filtering pass, and every new data set starts over. In exchange it
#          is a consistent estimator of the actual likelihood, with no
#          approximation beyond Monte Carlo error.
#
#   npe()  trains one conditional density estimator on simulations from the
#          prior. Training is the whole cost; after that, conditioning on any
#          data set of the same shape is a forward pass. In exchange the
#          posterior is an approximation whose error you have to check
#          (07_npe_diagnostics.R, 09_sbc.R).
#
# Model and data source
#   King, A. A., Ionides, E. L. and Breto, C. "Simulation-based Inference for
#   Epidemiological Dynamics" (SBIED), lesson 2: "Simulation of stochastic
#   dynamic models".
#   https://kingaa.github.io/sbied/stochsim/
#   R script: https://kingaa.github.io/sbied/stochsim/main.R
#   Data: https://kingaa.github.io/sbied/stochsim/Measles_Consett_1948.csv
#
#   The 1948 measles outbreak in Consett, England, population 38,000: 53 weeks
#   of case reports. The lesson's SIR model, in Euler steps of dt = 1/7 week,
#   with H accumulating recoveries within the week:
#     dN_SI ~ Binomial(S, 1 - exp(-Beta * I / N * dt))
#     dN_IR ~ Binomial(I, 1 - exp(-mu_IR * dt))
#     reports ~ NegBinomial(mu = rho * H, size = k)
#   initialised at S = round(eta * N), I = 1, R = round((1 - eta) * N).
#   The lesson simulates at Beta = 7.5, mu_IR = 0.5, rho = 0.5, k = 10,
#   eta = 0.03.
#
# The pomp half of this script is skipped if pomp is not installed. Everything
# above it runs regardless.
#
# Runtime: about 2 minutes, plus a few seconds the first time pomp compiles the
# model's C snippets.

library(neuralsbi)

has_torch <- requireNamespace("torch", quietly = TRUE) &&
  torch::torch_is_installed()
has_pomp <- requireNamespace("pomp", quietly = TRUE)

# ---------------------------------------------------------------------------
# The data
# ---------------------------------------------------------------------------

reports <- c(0, 0, 2, 0, 3, 0, 1, 0, 2, 4, 2, 4, 7, 34, 35, 22, 18, 75, 43,
             47, 44, 63, 49, 17, 19, 16, 1, 2, 0, 1, 1, 1, 1, 1, 4, 1, 0, 1,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
meas <- data.frame(week = seq_along(reports), reports = reports)
n_weeks <- nrow(meas)
N <- 38000
k_disp <- 10

# ---------------------------------------------------------------------------
# The simulator, in plain R
# ---------------------------------------------------------------------------
#
# This is the lesson's sir_step and sir_rinit written as one function of the
# parameters, returning one 53-week outbreak.

epidemic <- function(Beta, mu_IR, rho, eta) {
  S <- round(N * eta)
  I <- 1L
  out <- numeric(n_weeks)
  for (w in seq_len(n_weeks)) {
    H <- 0
    for (s in seq_len(7L)) {
      dN_SI <- rbinom(1L, S, 1 - exp(-Beta * I / N / 7))
      dN_IR <- rbinom(1L, I, 1 - exp(-mu_IR / 7))
      S <- S - dN_SI
      I <- I + dN_SI - dN_IR
      H <- H + dN_IR
    }
    out[w] <- rnbinom(1L, size = k_disp, mu = rho * H + 1e-6)
  }
  out
}

# As in 06_embedding_networks.R, the estimator sees log(1 + reports), because
# count noise scales with its mean and a single per-column standard deviation
# cannot cover a prior this wide. The particle filter below works on the raw
# counts, as it should, so the two methods are not handed the same numbers.
# They are handed the same model, which is the comparison that matters, and
# the log-likelihoods it reports are therefore directly comparable.
simulator <- function(Beta, mu_IR, rho, eta) {
  stats::setNames(log1p(epidemic(Beta, mu_IR, rho, eta)),
                  paste0("wk", seq_len(n_weeks)))
}

x_obs <- stats::setNames(log1p(reports), paste0("wk", seq_len(n_weeks)))

prior <- prior_uniform(
  low  = c(Beta =  1, mu_IR = 0.2, rho = 0.05, eta = 0.01),
  high = c(Beta = 60, mu_IR = 3.0, rho = 0.90, eta = 0.20)
)

# ---------------------------------------------------------------------------
# NPE
# ---------------------------------------------------------------------------
#
# An 8-feature embedding, because 53 weekly counts is a time series rather than
# a set of summaries. See 06_embedding_networks.R.

t0 <- Sys.time()
fit <- npe(prior, simulator, n_simulations = 4000,
           density_estimator = if (has_torch) "maf" else "linear_gaussian",
           embedding_net = if (has_torch) embedding_mlp(8L, c(64L, 64L)) else
             NULL,
           seed = 2024, verbose = TRUE)
npe_time <- as.numeric(Sys.time() - t0, units = "secs")
cat(sprintf("\nNPE: %.0f s for 4000 simulations plus training\n", npe_time))

post <- posterior(fit, x_obs = x_obs)
draws <- sample(post, 4000)
print(summary(draws))

npe_map <- map_estimate(post)
cat("\nMAP:\n"); print(round(npe_map, 3))
cat("SBIED lesson's simulation values: Beta 7.5, mu_IR 0.5, rho 0.5, eta 0.03\n")

# Amortization, priced. Conditioning on a second, third, hundredth outbreak
# costs a forward pass. pomp would need a full filtering run for each.
t0 <- Sys.time()
for (i in 1:20) invisible(sample(post, 1000, obs = x_obs))
cat(sprintf("20 further conditionings: %.2f s total\n",
            as.numeric(Sys.time() - t0, units = "secs")))

# ---------------------------------------------------------------------------
# pomp
# ---------------------------------------------------------------------------

if (!has_pomp) {
  cat("\npomp is not installed; skipping the comparison.\n")
  cat("  install.packages('pomp')\n")
} else {
  library(pomp)

  # The lesson's model, verbatim.
  sir_step <- Csnippet("
    double dN_SI = rbinom(S, 1 - exp(-Beta * I / N * dt));
    double dN_IR = rbinom(I, 1 - exp(-mu_IR * dt));
    S -= dN_SI;
    I += dN_SI - dN_IR;
    R += dN_IR;
    H += dN_IR;
  ")
  sir_rinit <- Csnippet("
    S = nearbyint(eta * N);
    I = 1;
    R = nearbyint((1 - eta) * N);
    H = 0;
  ")
  sir_dmeas <- Csnippet("lik = dnbinom_mu(reports, k, rho * H, give_log);")
  sir_rmeas <- Csnippet("reports = rnbinom_mu(k, rho * H);")

  measSIR <- pomp(
    meas, times = "week", t0 = 0,
    rprocess = euler(sir_step, delta.t = 1 / 7),
    rinit = sir_rinit, rmeasure = sir_rmeas, dmeasure = sir_dmeas,
    accumvars = "H",
    statenames = c("S", "I", "R", "H"),
    paramnames = c("Beta", "mu_IR", "N", "eta", "rho", "k")
  )

  pars <- function(theta) {
    c(Beta = unname(theta[["Beta"]]), mu_IR = unname(theta[["mu_IR"]]),
      rho = unname(theta[["rho"]]), eta = unname(theta[["eta"]]),
      k = k_disp, N = N)
  }

  # A particle-filter estimate of the log-likelihood, replicated so the Monte
  # Carlo error is visible. This is the number NPE never computes. The standard
  # error comes back NA when the replicates disagree so badly that one of them
  # dominates the log-mean-exp, which is itself the signal that the filter is
  # struggling at that parameter value.
  loglik_at <- function(theta, Np = 2000, reps = 5) {
    ll <- replicate(reps, logLik(pfilter(measSIR, params = pars(theta),
                                         Np = Np)))
    c(mean = logmeanexp(ll), se = logmeanexp(ll, se = TRUE)[["se"]])
  }

  t0 <- Sys.time()
  ll_lesson <- loglik_at(c(Beta = 7.5, mu_IR = 0.5, rho = 0.5, eta = 0.03))
  one_eval <- as.numeric(Sys.time() - t0, units = "secs")
  cat(sprintf("\none likelihood evaluation (5 x 2000 particles): %.1f s\n",
              one_eval))

  ll_map  <- loglik_at(npe_map)
  ll_mean <- loglik_at(colMeans(draws))

  cat("\nparticle-filter log-likelihood:\n")
  print(round(rbind(
    `SBIED simulation values` = ll_lesson,
    `NPE posterior mean`      = ll_mean,
    `NPE MAP`                 = ll_map), 2))

  # The lesson's values are a hand-chosen starting point for iterated
  # filtering, not a fit, so NPE's estimates landing above them is the outcome
  # to expect. If NPE lands well below, that is a real finding and the
  # diagnostics in 07 and 09 are where to take it.

  # How the pomp likelihood behaves across the NPE posterior. This is the
  # cross-check that costs the most and says the most: if the NPE posterior
  # concentrates where the true likelihood is high, the two methods agree
  # about the model even though only one of them ever wrote the likelihood
  # down.
  set.seed(1)
  idx <- base::sample.int(nrow(draws), 12)
  post_ll <- t(vapply(idx, function(i) loglik_at(draws[i, ], Np = 1000,
                                                 reps = 3), numeric(2)))
  spread <- data.frame(round(draws[idx, ], 3), loglik = round(post_ll[, 1], 1))
  print(spread[order(-spread$loglik), ])

  cat(sprintf(
    "\nbest of 12 posterior draws: %.1f;  NPE MAP: %.1f;  lesson values: %.1f\n",
    max(post_ll[, 1]), ll_map[["mean"]], ll_lesson[["mean"]]))

  # A short iterated-filtering run, for the local maximum-likelihood estimate.
  # Nmif = 20 is small; the SBIED lessons on mif2 use far more, from many
  # starting points. Treat this as a direction, not a converged answer.
  t0 <- Sys.time()
  mf <- mif2(measSIR,
             params = pars(npe_map),
             Np = 1000, Nmif = 20,
             cooling.fraction.50 = 0.5,
             rw.sd = rw_sd(Beta = 0.02, mu_IR = 0.02, rho = 0.02,
                           eta = ivp(0.02)),
             partrans = parameter_trans(log = c("Beta", "mu_IR"),
                                        logit = c("rho", "eta")),
             paramnames = c("Beta", "mu_IR", "rho", "eta", "N", "k"))
  cat(sprintf("\nmif2 (20 iterations, 1000 particles): %.0f s\n",
              as.numeric(Sys.time() - t0, units = "secs")))
  mif_par <- coef(mf)[c("Beta", "mu_IR", "rho", "eta")]
  print(round(mif_par, 3))
  print(round(loglik_at(mif_par), 2))

  cat("\ncomparison:\n")
  print(round(rbind(
    `NPE posterior mean` = colMeans(draws),
    `NPE MAP`            = npe_map,
    `mif2 from NPE MAP`  = mif_par), 3))
  print(round(rbind(
    `NPE posterior 2.5%`  = apply(draws, 2, quantile, 0.025),
    `NPE posterior 97.5%` = apply(draws, 2, quantile, 0.975)), 3))

  # mif2 returns a point, and its uncertainty needs a profile likelihood over
  # a grid: many more filtering runs. NPE returns the interval as part of the
  # same object, at no extra cost. That is the trade in one sentence.
}

# ---------------------------------------------------------------------------
# Where each one belongs
# ---------------------------------------------------------------------------
#
# pomp, when there is one data set, the likelihood is what you want, and you
# are willing to pay per evaluation for it. Iterated filtering with a profile
# likelihood is the mature, well-understood answer for a single POMP.
#
# npe(), when the same model has to be fit to many data sets (51 jurisdictions,
# 200 hospitals, a simulation study), or when you want a posterior rather than
# a point estimate, or when the simulator is a black box with no filtering
# structure to exploit at all. The training cost is paid once.
#
# Both, when it matters. A particle filter evaluated at the NPE posterior mode
# is the cheapest independent check on an NPE fit that exists for a POMP, and
# it does not depend on any assumption the estimator made.
