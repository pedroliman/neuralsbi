# Posterior log-density

Posterior log-density

## Usage

``` r
log_prob(post, theta, x = NULL, normalize = TRUE, n_normalization = 10000L)
```

## Arguments

- post:

  An `nsbi_posterior` object.

- theta:

  Matrix (or vector) of parameter values to evaluate.

- x:

  Observation to condition on (defaults to `x_obs`).

- normalize:

  For bounded priors, renormalize by the estimated acceptance
  probability and return `-Inf` outside the prior support.

- n_normalization:

  Number of draws used to estimate the normalizing (acceptance) constant
  when `normalize = TRUE`.

## Value

Numeric vector of log posterior densities.
