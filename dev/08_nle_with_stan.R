# 08_nle_with_stan.R ---------------------------------------------------------
#
# A learned likelihood is only half useful while it is trapped inside R.
# stan_code() transpiles an nle() fit into a Stan functions block: the network
# weights travel as data, Stan differentiates the generated code itself, and
# nothing has to link against torch at run time. The point is not to replace
# the package's own slice sampler. It is that the surrogate stops being the
# whole model and becomes one term in a model you write, with hierarchy,
# covariates, or a second data source whose likelihood you do know.
#
# Model and data source
#   Viechtbauer, W. metafor / metadat. Dataset dat.bcg (now
#   metadat::dat.colditz1994), the 13 randomized trials of the BCG vaccine
#   against tuberculosis analysed in
#   Colditz, G. A. et al. (1994), "Efficacy of BCG vaccine in the prevention of
#   tuberculosis", JAMA 271(9), 698-702.
#   https://wviechtb.github.io/metadat/reference/dat.colditz1994.html
#
#   metafor's own worked example fits the random-effects model
#     dat <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
#                   data = dat.bcg)
#     rma(yi, vi, data = dat)
#   https://wviechtb.github.io/metafor/reference/rma.uni.html
#
#   The log risk ratios and their standard errors below are computed from that
#   dataset with the usual formulas. Direct maximum likelihood on the exact
#   random-effects likelihood gives mu = -0.711, tau = 0.529, which is the
#   number to hold the surrogate against.
#
# Runtime: about 2 minutes, plus a Stan compile if cmdstanr or rstan is
# installed.

library(neuralsbi)

has_torch <- requireNamespace("torch", quietly = TRUE) &&
  torch::torch_is_installed()

# ---------------------------------------------------------------------------
# The data
# ---------------------------------------------------------------------------

yi <- c(-0.8893, -1.5854, -1.3481, -1.4416, -0.2175, -0.7861, -1.6209,
         0.0120, -0.4694, -1.3713, -0.3394,  0.4459, -0.0173)
si <- c( 0.5706,  0.4411,  0.6445,  0.1415,  0.2263,  0.0831,  0.4722,
         0.0629,  0.2376,  0.2702,  0.1114,  0.7297,  0.2672)

# One observation is one trial. The within-trial standard error is part of what
# is observed, so it goes into x rather than being pretended away: each row is
# (log risk ratio, log standard error). Since s does not depend on the
# parameters, its contribution to the likelihood is a constant that cancels in
# the posterior, and the trials become exchangeable, which is what nle() needs
# to sum over them.
x_obs <- cbind(y = yi, log_s = log(si))
print(round(x_obs, 3))

# ---------------------------------------------------------------------------
# The model
# ---------------------------------------------------------------------------
#
#   s_i ~ empirical distribution of the observed standard errors
#   y_i | s_i ~ Normal(mu, sqrt(tau^2 + s_i^2))

simulator <- function(mu, log_tau) {
  s <- si[base::sample.int(length(si), 1L)]
  c(y = rnorm(1, mu, sqrt(exp(2 * log_tau) + s^2)), log_s = log(s))
}

prior <- prior_uniform(low  = c(mu = -2.0, log_tau = log(0.02)),
                       high = c(mu =  1.0, log_tau = log(2.00)))

# "mdn" and "maf" and "linear_gaussian" export to Stan. "nsf" does not: its
# rational-quadratic spline would be a large and fragile block of generated
# code. Refit with "maf" if you need a flow.
estimator <- if (has_torch) "mdn" else "linear_gaussian"

fit <- nle(prior, simulator, n_simulations = 5000,
           density_estimator = estimator, seed = 2024, verbose = TRUE)
print(fit)

# ---------------------------------------------------------------------------
# Does the surrogate know the likelihood?
# ---------------------------------------------------------------------------
#
# This model has a closed form, which makes it a good place to check the
# machinery. Compare the surrogate's log-likelihood surface against the exact
# one over a grid of mu.

exact_loglik <- function(mu, log_tau) {
  sum(dnorm(yi, mu, sqrt(exp(2 * log_tau) + si^2), log = TRUE))
}

mu_grid <- seq(-1.4, 0.0, length.out = 9)
lt <- log(0.53)
comp <- data.frame(
  mu = mu_grid,
  exact = vapply(mu_grid, exact_loglik, numeric(1), log_tau = lt),
  surrogate = log_lik(fit, cbind(mu = mu_grid, log_tau = lt), x_obs)
)
# Both are defined up to a constant that does not depend on mu (the surrogate
# also models log_s), so compare after centring on each column's maximum.
comp$exact <- comp$exact - max(comp$exact)
comp$surrogate <- comp$surrogate - max(comp$surrogate)
print(round(comp, 2))

cat("\ncorrelation of the two surfaces:",
    sprintf("%.4f", cor(comp$exact, comp$surrogate)), "\n")

# Expect the shape to match and individual points to wobble by a nat or two.
# That is what a mixture density network trained on five thousand simulations
# buys: the right curvature, not the exact number. Raise n_simulations if the
# wobble matters, and use it as the check when there is no closed form to
# compare against, by holding out simulations instead.

# ---------------------------------------------------------------------------
# The posterior, before Stan gets involved
# ---------------------------------------------------------------------------

post <- posterior(fit, x_obs, n_chains = 20, warmup = 200, thin = 2, seed = 5)
draws <- sample(post, 3000)
print(summary(draws))
cat("\nexact-likelihood MLE: mu = -0.711, tau = 0.529 (log_tau = -0.636)\n")

# ---------------------------------------------------------------------------
# stan_code(): the generated program
# ---------------------------------------------------------------------------
#
# Two entry points come out, where `w` is the packed weight vector:
#
#   real nsbi_log_lik_lpdf(vector x, vector theta, vector w);      // one row
#   real nsbi_log_lik_sum_lpdf(matrix x, vector theta, vector w);  // i.i.d. rows
#
# Both take x and theta in the original units. The standardization the
# estimator trained under is folded into the generated code, and so is its
# Jacobian.

code <- stan_code(fit)
lines <- strsplit(code, "\n")[[1]]
cat("\ngenerated program:", length(lines), "lines\n\n")
cat(paste(head(lines, 40), collapse = "\n"), "\n  ...\n")

# The model blocks at the end restate the prior as a Stan sampling statement,
# which works because prior_uniform() and prior_normal() are named
# distributions with parameters. A prior_custom() is arbitrary R code and
# cannot be restated, so use model = FALSE and write the model block yourself.
data_at <- grep("^data \\{", lines)[1]
cat("\n", paste(lines[data_at:length(lines)], collapse = "\n"), "\n")

# functions only, for #include-ing into a model of your own
fns <- stan_code(fit, model = FALSE)
cat("\nfunctions block alone:", length(strsplit(fns, "\n")[[1]]), "lines\n")

# ---------------------------------------------------------------------------
# stan_data(): the weights and the observation, ready to pass
# ---------------------------------------------------------------------------

sdata <- stan_data(fit, x_obs)
str(sdata, max.level = 1)

# nsbi_w is the packed weight vector, nsbi_nw its length, N and x the
# observation, and nsbi_low/nsbi_high the uniform prior's box.

# write_stan_model() puts the program on disk.
stan_file <- file.path(tempdir(), "bcg_surrogate.stan")
write_stan_model(fit, stan_file)
cat("wrote", stan_file, "\n")

# ---------------------------------------------------------------------------
# Running it
# ---------------------------------------------------------------------------
#
# posterior(fit, sampler = "stan") does the whole round trip: generate,
# compile, run NUTS, return draws. It needs cmdstanr or rstan, and pays a
# one-time compile. NUTS mixes better than the slice sampler on correlated
# posteriors, which is the reason to bother.

# Check for CmdStan itself, not just the cmdstanr R package.
# install.packages("cmdstanr") does not install CmdStan, and a session with the
# R package but no toolchain fails inside cmdstan_model() with "CmdStan path
# has not been set".
has_cmdstan <- requireNamespace("cmdstanr", quietly = TRUE) &&
  !is.null(tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE),
                    error = function(e) NULL))
has_rstan <- requireNamespace("rstan", quietly = TRUE)

if (has_cmdstan || has_rstan) {
  post_stan <- posterior(fit, x_obs, sampler = "stan",
                         n_chains = 4, iter_warmup = 500, iter_sampling = 500,
                         seed = 5)
  draws_stan <- sample(post_stan, 2000)
  print(summary(draws_stan))
  # The two samplers target the same distribution, so their draws should be
  # hard to tell apart.
  print(c2st(draws, draws_stan, seed = 1))
} else {
  cat("\nNo working Stan back end: skipping the NUTS run.\n")
  cat("  install.packages('cmdstanr',",
      "repos = c('https://stan-dev.r-universe.dev', getOption('repos')))\n")
  cat("  cmdstanr::install_cmdstan()\n")
}

# ---------------------------------------------------------------------------
# The reason to generate code rather than call a callable
# ---------------------------------------------------------------------------
#
# Here is the model the export is for. The surrogate supplies the likelihood of
# one trial's result. Everything else is a model you write: a hierarchical
# prior on mu across subgroups, a covariate, a second evidence source. None of
# that is reachable from a posterior estimator, which only knows the one
# conditional it was trained on.

sketch <- '
functions {
  #include bcg_surrogate_functions.stan   // from stan_code(fit, model = FALSE)
}
data {
  int<lower=1> nsbi_nw;   vector[nsbi_nw] nsbi_w;   // trained weights
  int<lower=1> N;         matrix[N, 2] x;           // trials: (log RR, log SE)
  int<lower=1> K;         array[N] int<lower=1, upper=K> region;
}
parameters {
  real mu_global;
  real<lower=0> sigma_region;
  vector[K] mu_region;
  real log_tau;
}
model {
  // a hierarchy over regions, which the surrogate knows nothing about
  mu_global    ~ normal(0, 1);
  sigma_region ~ exponential(2);
  mu_region    ~ normal(mu_global, sigma_region);
  log_tau      ~ uniform(log(0.02), log(2));

  // the learned likelihood, one trial at a time, each under its own region mean
  for (n in 1:N)
    target += nsbi_log_lik_lpdf(to_vector(x[n]) |
                                [mu_region[region[n]], log_tau]\', nsbi_w);
}
'
cat(sketch)

# ---------------------------------------------------------------------------
# Limits worth knowing
# ---------------------------------------------------------------------------
#
#   * Only nle() fits export. An nre() fit holds a classifier, not a density,
#     and there is no p(x | theta) in it to write out.
#   * "nsf" does not export. Use "maf" or "mdn".
#   * A fit restored with readRDS() has a dead network pointer, and stan_code()
#     says so rather than failing later. Use save_nle()/load_nle().
#   * posterior(sampler = "stan") prefers cmdstanr whenever the R package is
#     installed. If CmdStan itself is not installed it fails there rather than
#     falling back to rstan, so run cmdstanr::install_cmdstan() first or
#     uninstall cmdstanr if rstan is what you have.
#   * The generated code is source you can read and edit. That is deliberate:
#     the package's own tests evaluate the emitted functions and compare them
#     against log_lik().
