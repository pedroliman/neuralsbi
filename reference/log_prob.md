# Posterior log-density

Posterior log-density

## Usage

``` r
# S3 method for class 'nsbi_nle_posterior'
log_prob(post, theta, x = NULL, normalize = TRUE, ...)

log_prob(post, theta, x = NULL, ...)

# S3 method for class 'nsbi_posterior'
log_prob(
  post,
  theta,
  x = NULL,
  normalize = TRUE,
  n_normalization = 10000L,
  ...
)
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

- ...:

  Passed to methods.

- n_normalization:

  Number of draws used to estimate the normalizing (acceptance) constant
  when `normalize = TRUE`.

## Value

Numeric vector of log posterior densities. For a posterior built from an
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) fit the
value is **unnormalized** – the evidence \\p(x)\\ is not available – so
differences between two `theta` are meaningful but the absolute level is
not.
