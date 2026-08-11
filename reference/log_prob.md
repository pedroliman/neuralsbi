# Posterior log-density

Posterior log-density

## Usage

``` r
# S3 method for class 'nsbi_mcmc_posterior'
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
  when `normalize = TRUE`. If none of them land inside the prior
  support, the estimate is floored at `1 / n_normalization` to avoid
  `log(0)` and a warning says so – the same warning
  [`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md)
  raises when rejection sampling comes up empty.

## Value

Numeric vector of log posterior densities. For a posterior built from an
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) or
[`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md) fit the
value is **unnormalized** – the evidence \\p(x)\\ is not available – so
differences between two `theta` are meaningful but the absolute level is
not.
