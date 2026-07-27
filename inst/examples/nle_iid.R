# Neural likelihood estimation with repeated independent observations.
#
# The model is the g-and-k distribution: defined by its quantile function, so
# simulating is one line and the density has no closed form. That combination
# is the standard test case for likelihood-free inference, and the repeated
# observations are what make NLE the right method rather than NPE.
#
# Run with:  Rscript inst/examples/nle_iid.R

library(neuralsbi)
set.seed(1)

## The simulator ------------------------------------------------------------

rgk <- function(u, A, B, g, k, c = 0.8) {
  z <- qnorm(u)
  A + B * (1 + c * (1 - exp(-g * z)) / (1 + exp(-g * z))) * (1 + z^2)^k * z
}

# One call, one observation.
simulator <- function(A, B, g, k) c(y = rgk(runif(1), A, B, g, k))

prior <- prior_uniform(low  = c(A = 0, B = 0, g = 0, k = 0),
                       high = c(A = 10, B = 10, g = 10, k = 10))

## Train once ---------------------------------------------------------------

fit <- nle(prior, simulator, n_simulations = 20000,
           density_estimator = "maf", seed = 1, verbose = TRUE)
print(fit)

## Condition on however many observations you have --------------------------

theta_true <- c(A = 3, B = 1, g = 2, k = 0.5)
x_all <- matrix(rgk(runif(5000), 3, 1, 2, 0.5), ncol = 1)

for (n_obs in c(50, 500, 5000)) {
  post <- posterior(fit, x_all[seq_len(n_obs), , drop = FALSE],
                    n_chains = 20, warmup = 200, thin = 10, seed = 2)
  draws <- sample(post, 4000)
  cat(sprintf("\nn = %d\n", n_obs))
  print(round(rbind(mean = colMeans(draws),
                    sd = apply(draws, 2, sd),
                    truth = theta_true), 3))
  cat(sprintf("max Rhat %.3f, min bulk ESS %.0f\n",
              max(attr(draws, "diagnostics")$rhat),
              min(attr(draws, "diagnostics")$ess_bulk)))
}

## The same likelihood, as a plain R function -------------------------------

loglik <- likelihood_fn(fit, x_all[1:500, , drop = FALSE])
cat("\nlog-likelihood at the truth:", loglik(theta_true), "\n")

## Or as Stan code ----------------------------------------------------------

path <- file.path(tempdir(), "gk_likelihood.stan")
write_stan_model(fit, path)
cat("Stan model written to", path, "\n")

if (requireNamespace("cmdstanr", quietly = TRUE)) {
  by_stan <- sample(
    posterior(fit, x_all[1:500, , drop = FALSE], sampler = "stan",
              n_chains = 4, seed = 3),
    4000)
  cat("\nNUTS posterior mean:", round(colMeans(by_stan), 3), "\n")
}
