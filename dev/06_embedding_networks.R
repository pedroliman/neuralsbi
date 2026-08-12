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
#   Measurement: reports ~ NegBinomial(mu = rho * H, size = k).
#   The lesson simulates at Beta = 7.5, mu_IR = 0.5, rho = 0.5, k = 10,
#   eta = 0.03.
#
# Runtime: about 4 minutes on a laptop CPU.

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

simulator <- function(Beta, mu_IR, rho, eta) {
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
  stats::setNames(out, paste0("wk", seq_len(n_weeks)))
}

prior <- prior_uniform(
  low  = c(Beta =  1, mu_IR = 0.2, rho = 0.05, eta = 0.005),
  high = c(Beta = 30, mu_IR = 3.0, rho = 0.90, eta = 0.100)
)

# ---------------------------------------------------------------------------
# Simulate once, train twice
# ---------------------------------------------------------------------------
#
# The comparison is only fair if both estimators see the same simulations.

sims <- simulate_for_sbi(simulator, prior, n = 3000, seed = 11)
cat("simulated", nrow(sims$x), "epidemics of", ncol(sims$x), "weeks;",
    sims$n_dropped, "dropped\n")

# ---------------------------------------------------------------------------
# 1. No embedding: the flow conditions on all 53 weeks
# ---------------------------------------------------------------------------

t0 <- Sys.time()
fit_raw <- npe(prior, theta = sims$theta, x = sims$x,
               density_estimator = "maf", seed = 1)
cat(sprintf("no embedding: %.1f s\n",
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
fit_emb <- npe(prior, theta = sims$theta, x = sims$x,
               density_estimator = "maf", embedding_net = emb, seed = 1)
cat(sprintf("with embedding: %.1f s\n",
            as.numeric(Sys.time() - t0, units = "secs")))
info_emb <- summary(fit_emb)

# The embedding sees the standardized data, the same z-scoring npe() applies to
# x without one, so its inputs are on a common scale. Standardizing the
# features is left to the network.

cat("\nbest validation loss (lower is better):\n")
cat(sprintf("  no embedding  : %.3f\n", info_raw$best_val_loss))
cat(sprintf("  with embedding: %.3f\n", info_emb$best_val_loss))

# ---------------------------------------------------------------------------
# Posteriors, side by side
# ---------------------------------------------------------------------------

post_raw <- posterior(fit_raw, x_obs = reports)
post_emb <- posterior(fit_emb, x_obs = reports)

d_raw <- sample(post_raw, 3000)
d_emb <- sample(post_emb, 3000)

cat("\nno embedding:\n");   print(summary(d_raw))
cat("\nwith embedding:\n"); print(summary(d_emb))

# c2st() puts a number on how different the two posteriors are. Near 0.5 means
# a linear classifier cannot tell them apart.
print(c2st(d_raw, d_emb, seed = 1))

# ---------------------------------------------------------------------------
# Which one reproduces the outbreak?
# ---------------------------------------------------------------------------
#
# The honest test is predictive, not internal to the estimator. Push posterior
# draws back through the simulator and see which set of predictions contains
# the observed series.

pp <- function(post, label) {
  pred <- posterior_predictive(post, simulator, n = 300)
  band <- apply(pred, 2, quantile, probs = c(0.05, 0.95))
  inside <- mean(reports >= band[1, ] & reports <= band[2, ])
  peak_obs <- max(reports)
  peak_pred <- apply(pred, 1, max)
  cat(sprintf("%-15s weeks inside 90%% band: %.0f%%   peak %d vs predicted %.0f [%.0f, %.0f]\n",
              label, 100 * inside, peak_obs, median(peak_pred),
              quantile(peak_pred, 0.05), quantile(peak_pred, 0.95)))
  invisible(pred)
}
cat("\n")
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

fit_tiny <- npe(prior, theta = sims$theta, x = sims$x,
                density_estimator = "maf",
                embedding_net = embedding_mlp(output_dim = 2L,
                                              hidden = c(32L, 32L)),
                seed = 1)
cat(sprintf("\noutput_dim = 2: val loss %.3f\n",
            summary(fit_tiny)$best_val_loss))
print(summary(sample(posterior(fit_tiny, x_obs = reports), 2000)))

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
