# Shared across test files: the conjugate-Gaussian reference posterior for
# theta ~ N(0, I), x | theta ~ N(theta, sigma^2 I). Every linear_gaussian
# parity test and every neural-estimator smoke test (MAF, MDN, NSF,
# embedding, sequential) checks its draws against this closed form, so it is
# the one reference the whole Level 1 verification strategy rests on.
analytic_gauss_posterior <- function(x_obs, sigma, d) {
  prec <- diag(d) + diag(d) / sigma^2
  Sigma <- solve(prec)
  mu <- as.numeric(Sigma %*% (x_obs / sigma^2))
  list(mu = mu, Sigma = Sigma)
}
