# 06_embedding_networks.R ----------------------------------------------------
#
# When the observation is a time series rather than a handful of summaries,
# feeding all of it straight into the density estimator spends the network's
# capacity on the wrong thing. An embedding network learns a low-dimensional
# summary f(x) jointly with the estimator, so the conditional becomes
# q(theta | f(x)). embedding_mlp() builds one; pass it to npe() or nre() via
# embedding_net.
#
# Model and data source
#   King, A. A., Ionides, E. L. and Breto, C. "Simulation-based Inference for
#   Epidemiological Dynamics" (SBIED), lesson 2: "Simulation of stochastic
#   dynamic models".
#   https://kingaa.github.io/sbied/stochsim/
#   R script: https://kingaa.github.io/sbied/stochsim/main.R
#
#   The lesson's SIR model for the 1948 measles outbreak in Consett, England
#   (population 38,000), fit with pomp. Weekly case reports are
#   Measles_Consett_1948.csv, reproduced below.
#   https://kingaa.github.io/sbied/stochsim/Measles_Consett_1948.csv
#
#   Process, in Euler steps of dt = 1/7 week:
#     dN_SI ~ Binomial(S, 1 - exp(-Beta * I / N * dt))
#     dN_IR ~ Binomial(I, 1 - exp(-mu_IR * dt))
#   with H accumulating dN_IR within each week and resetting at the week
#   boundary. Initial state S = round(eta * N), I = 1, R = round((1 - eta) * N).
#   Measurement: reports ~ NegBinomial(mu = rho * H, size = k). k is fixed at
#   the lesson's 10 and the other four parameters are estimated.
#
# Runtime: about 5 minutes on a laptop CPU.

library(neuralsbi)

if (!(requireNamespace("torch", quietly = TRUE) &&
      torch::torch_is_installed())) {
  stop("This script needs torch: embedding networks are a neural feature.")
}

# ---------------------------------------------------------------------------
# The data
# ---------------------------------------------------------------------------

reports <- c(0, 0, 2, 0, 3, 0, 1, 0, 2, 4, 2, 4, 7, 34, 35, 22, 18, 75, 43,
             47, 44, 63, 49, 17, 19, 16, 1, 2, 0, 1, 1, 1, 1, 1, 4, 1, 0, 1,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
n_weeks <- length(reports)
N <- 38000
k_disp <- 10        # measurement overdispersion, fixed at the lesson's value

# ---------------------------------------------------------------------------
# The simulator
# ---------------------------------------------------------------------------

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

# The simulator returns log(1 + reports). Weekly counts here run from 0 to
# several hundred across the prior, and z-scoring a series that skewed leaves
# the estimator fitting the few large weeks and ignoring the shape. Fitting on
# the log scale is worth more here than any change of architecture, which is
# the general lesson: try the transform before the bigger network.
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
# What this model does before any inference happens
# ---------------------------------------------------------------------------
#
# The chain starts at I = 1, so a large share of prior draws produce no
# outbreak at all: the index case recovers before infecting anyone. Knowing
# that fraction changes how you read everything below, so measure it first.

sims <- simulate_for_sbi(simulator, prior, n = 4000, seed = 11)
peaks <- apply(expm1(sims$x), 1, max)
cat(sprintf("prior predictive: %.0f%% of draws never exceed 2 cases in a week\n",
            100 * mean(peaks <= 2)))
cat(sprintf("observed peak: %d cases in week %d, %d cases in total\n",
            max(reports), which.max(reports), sum(reports)))

# ---------------------------------------------------------------------------
# 1. No embedding: the flow conditions on all 53 weeks
# ---------------------------------------------------------------------------
#
# max_epochs and patience are trimmed so the script finishes in a few minutes.
# All three fits get the same deal and the same simulations, which is what
# makes the comparison about the embedding.

ctl <- list(max_epochs = 200L, patience = 12L)

t0 <- Sys.time()
fit_raw <- do.call(npe, c(list(prior, theta = sims$theta, x = sims$x,
                               density_estimator = "maf", seed = 1), ctl))
cat(sprintf("\nno embedding: %.0f s\n",
            as.numeric(Sys.time() - t0, units = "secs")))
info_raw <- summary(fit_raw)

# ---------------------------------------------------------------------------
# 2. With an embedding: 53 weeks compressed to 8 learned features
# ---------------------------------------------------------------------------
#
# embedding_mlp() is a specification, not a network. It carries no torch
# objects, so it can be built without torch installed; the module itself is
# constructed at fit time and its weights live inside the fitted estimator.
# Sampling and log_prob route through it automatically.

emb <- embedding_mlp(output_dim = 8L, hidden = c(64L, 64L))
str(emb)

t0 <- Sys.time()
fit_emb <- do.call(npe, c(list(prior, theta = sims$theta, x = sims$x,
                               density_estimator = "maf",
                               embedding_net = emb, seed = 1), ctl))
cat(sprintf("with embedding: %.0f s\n",
            as.numeric(Sys.time() - t0, units = "secs")))
info_emb <- summary(fit_emb)

# The embedding sees the standardized data, the same z-scoring npe() applies to
# x without one, so its inputs are on a common scale. Standardizing the
# features is left to the network.

cat("\nbest validation loss (lower is better):\n")
cat(sprintf("  no embedding  : %.3f  (%d epochs)\n",
            info_raw$best_val_loss, info_raw$epochs_trained))
cat(sprintf("  with embedding: %.3f  (%d epochs)\n",
            info_emb$best_val_loss, info_emb$epochs_trained))

# ---------------------------------------------------------------------------
# Posteriors, side by side
# ---------------------------------------------------------------------------

post_raw <- posterior(fit_raw, x_obs = x_obs)
post_emb <- posterior(fit_emb, x_obs = x_obs)

d_raw <- sample(post_raw, 3000)
d_emb <- sample(post_emb, 3000)

cat("\nno embedding:\n");   print(summary(d_raw))
cat("\nwith embedding:\n"); print(summary(d_emb))

# c2st() puts a number on how different the two posteriors are. Near 0.5 means
# a linear classifier cannot tell them apart.
print(c2st(d_raw, d_emb, seed = 1))

# The four parameters are not separately identified by one outbreak: Beta, eta
# and mu_IR enter the early growth almost only through Beta * eta / mu_IR, and
# the prior stays visible in each of them on its own. The combination is what
# the data pin down, so that is what to compare.
r_eff <- function(d) d[, "Beta"] * d[, "eta"] / d[, "mu_IR"]
cat(sprintf("\ninitial R_eff = Beta * eta / mu_IR\n"))
cat(sprintf("  prior         : %.2f [%.2f, %.2f]\n",
            median(r_eff(sims$theta)), quantile(r_eff(sims$theta), 0.05),
            quantile(r_eff(sims$theta), 0.95)))
cat(sprintf("  no embedding  : %.2f [%.2f, %.2f]\n",
            median(r_eff(d_raw)), quantile(r_eff(d_raw), 0.05),
            quantile(r_eff(d_raw), 0.95)))
cat(sprintf("  with embedding: %.2f [%.2f, %.2f]\n",
            median(r_eff(d_emb)), quantile(r_eff(d_emb), 0.05),
            quantile(r_eff(d_emb), 0.95)))

# ---------------------------------------------------------------------------
# Which one reproduces the outbreak?
# ---------------------------------------------------------------------------
#
# Read this one carefully. The model has a stochastic extinction mode at I = 1,
# so even at the right parameters most predictive draws fizzle out. A raw
# predictive band therefore looks terrible for any fit, and would be just as
# terrible if we knew the parameters exactly. What the data condition on is an
# outbreak having happened, so condition the check the same way and compare the
# outbreaks with the outbreak.

pp <- function(post, label) {
  pred <- expm1(posterior_predictive(post, simulator, n = 400))
  took_off <- apply(pred, 1, max) > 10
  cat(sprintf("%-15s %.0f%% of draws take off; among those, peak %.0f [%.0f, %.0f] and total %.0f [%.0f, %.0f]\n",
              label, 100 * mean(took_off),
              median(apply(pred[took_off, ], 1, max)),
              quantile(apply(pred[took_off, ], 1, max), 0.05),
              quantile(apply(pred[took_off, ], 1, max), 0.95),
              median(rowSums(pred[took_off, ])),
              quantile(rowSums(pred[took_off, ]), 0.05),
              quantile(rowSums(pred[took_off, ]), 0.95)))
  invisible(pred)
}
cat(sprintf("\nobserved: peak %d, total %d\n", max(reports), sum(reports)))
pp(post_raw, "no embedding")
pp(post_emb, "with embedding")

# ---------------------------------------------------------------------------
# Choosing output_dim
# ---------------------------------------------------------------------------
#
# output_dim is the effective data dimension the estimator conditions on. Too
# small and the summary throws away information the parameters need; too large
# and you are back to the original problem with extra layers in front of it.
# With four parameters, an 8-feature summary is already generous. Here is what
# a very tight one does.

fit_tiny <- do.call(npe, c(list(prior, theta = sims$theta, x = sims$x,
                                density_estimator = "maf",
                                embedding_net = embedding_mlp(2L, c(32L, 32L)),
                                seed = 1), ctl))
info_tiny <- summary(fit_tiny)
cat(sprintf("\noutput_dim = 2: val loss %.3f (%d epochs)\n",
            info_tiny$best_val_loss, info_tiny$epochs_trained))
d_tiny <- sample(posterior(fit_tiny, x_obs = x_obs), 3000)
cat(sprintf("  R_eff %.2f [%.2f, %.2f]\n", median(r_eff(d_tiny)),
            quantile(r_eff(d_tiny), 0.05), quantile(r_eff(d_tiny), 0.95)))

# ---------------------------------------------------------------------------
# Notes
# ---------------------------------------------------------------------------
#
# nre() takes embedding_net too, and the classifier then sees (theta, f(x)).
#
# The "linear_gaussian" estimator and the "logistic" classifier ignore
# embedding_net and warn, because neither has a network to train it inside.
#
# An embedding is not a substitute for good summaries when you have them. If
# you already know that the peak height and the epidemic duration are what
# identify the parameters, passing those two numbers as x is cheaper and more
# transparent than learning them.
#
# And it is not a substitute for a model whose parameters are identified. No
# summary network can separate Beta from eta here, because the data do not.
# 15_npe_vs_pomp.R takes the same model to a particle filter, which agrees.
