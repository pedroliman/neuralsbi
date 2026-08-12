# 00_installation.R ----------------------------------------------------------
#
# Check that neuralsbi and its optional back ends are working before you spend
# a simulation budget on them. Nothing here fits a real model: it is the
# five-minute check you run once, then move on to 01_basic_npe_example.R.
#
# No external simulator source: this script has no scientific content.
#
# Runtime: a few seconds (or a few minutes the first time, if libtorch has to
# download).

# ---------------------------------------------------------------------------
# 1. Install the package
# ---------------------------------------------------------------------------

# From CRAN, once released:
#   install.packages("neuralsbi")
#
# Development version:
#   # install.packages("pak")
#   pak::pak("pedroliman/neuralsbi")
#
# From a local clone (what the rest of these scripts assume):
#   devtools::load_all()   # interactive work
#   R CMD INSTALL .        # or a real install

library(neuralsbi)

cat("neuralsbi", as.character(packageVersion("neuralsbi")), "\n")
cat("R", R.version.string, "\n\n")

# ---------------------------------------------------------------------------
# 2. The torch back end
# ---------------------------------------------------------------------------
#
# torch is a Suggests, not an Imports. The package loads and the
# "linear_gaussian" estimator works without it. The three neural estimators
# ("maf", "mdn", "nsf"), embedding networks and neural ratio estimation all
# need it.
#
# Installing torch is two steps: the R package, then the libtorch C++ library
# it binds to.
#
#   install.packages("torch")
#   torch::install_torch()      # downloads libtorch, ~200 MB

has_torch <- requireNamespace("torch", quietly = TRUE) &&
  torch::torch_is_installed()

if (has_torch) {
  cat("torch:", as.character(packageVersion("torch")),
      "with libtorch installed\n")
  # A forward pass, to be sure the shared library actually loads.
  z <- torch::torch_randn(2, 2)
  cat("torch smoke test:", is.numeric(as.matrix(z)), "\n\n")
} else {
  cat("torch is NOT available.\n")
  cat("  install.packages('torch'); torch::install_torch()\n")
  cat("  Until then, use density_estimator = 'linear_gaussian'.\n\n")
}

# ---------------------------------------------------------------------------
# 3. The other optional dependencies
# ---------------------------------------------------------------------------
#
# Each one buys a specific feature. None of them is needed to fit a model.

optional <- c(
  ggplot2   = "all plot_*() functions",
  GGally    = "pairplot()",
  ggdensity = "pairplot() density regions",
  future    = "parallel simulation (see 14_utilities.R)",
  progressr = "progress bars with an ETA",
  cmdstanr  = "posterior(fit, sampler = 'stan') and running exported models",
  rstan     = "the same, via rstan instead of cmdstanr",
  posterior = "MCMC diagnostics reported by the Stan sampler"
)

status <- vapply(names(optional), requireNamespace, logical(1), quietly = TRUE)
print(data.frame(package = names(optional),
                 installed = status,
                 buys = unname(optional),
                 row.names = NULL))
cat("\n")

# ---------------------------------------------------------------------------
# 4. Smoke test: a two-parameter fit with no torch involved
# ---------------------------------------------------------------------------
#
# "linear_gaussian" is a closed-form conditional Gaussian. It is exact when the
# model really is linear-Gaussian, which makes it both the fastest way to check
# that the pipeline runs end to end and the regression oracle the package tests
# itself against.

prior <- prior_uniform(low  = c(mu = -3, nu = -3),
                       high = c(mu =  3, nu =  3))

simulator <- function(mu, nu) {
  c(a = mu + rnorm(1, sd = 0.2),
    b = nu + rnorm(1, sd = 0.2))
}

fit <- npe(prior, simulator, n_simulations = 1000,
           density_estimator = "linear_gaussian", seed = 1)
print(fit)

post <- posterior(fit, x_obs = c(a = 1.0, b = -0.5))
draws <- sample(post, 2000)
print(summary(draws))

# The posterior mean should sit near the observation, since the simulator is
# the identity plus small noise and the prior is flat over a wide box.
cat("\nposterior means:", sprintf("%.3f", colMeans(draws)), "\n")
cat("observation    :", c(1.0, -0.5), "\n")

# ---------------------------------------------------------------------------
# 5. One thing to know about sample()
# ---------------------------------------------------------------------------
#
# neuralsbi defines sample() as an S3 generic, which masks base::sample. The
# default method forwards to base::sample, so ordinary use is unaffected:

set.seed(1)
print(sample(1:10, 3))

# But inside a package or a script where the masking matters, call
# base::sample() explicitly, or use the non-generic alias for posterior draws:
draws2 <- sample_posterior(post, n = 100)
cat("sample_posterior() gave", nrow(draws2), "draws\n")

# ---------------------------------------------------------------------------
# Where to go next
# ---------------------------------------------------------------------------
#
#   01_basic_npe_example.R   the shortest end-to-end fit on real data
#   02_priors.R              building priors
#   03_npe.R                 every NPE option worth knowing
#   11_utilities.R and up    saving fits, parallel simulation, diagnostics
