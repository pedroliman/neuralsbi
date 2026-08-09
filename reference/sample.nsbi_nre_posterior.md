# Sample an NRE posterior with MCMC

Identical in behavior to
[`sample.nsbi_nle_posterior()`](https://neuralsbi.pedrodelima.com/reference/sample.nsbi_nle_posterior.md),
including the draw cache: asking for the same or fewer draws from the
same observation returns the cached run instead of starting a new chain.

## Usage

``` r
# S3 method for class 'nsbi_nre_posterior'
sample(
  x,
  size = 1000,
  n = size,
  obs = NULL,
  refresh = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  An `nsbi_nre_posterior` object (named `x` to satisfy the
  [`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md)
  generic).

- size, n:

  Number of posterior draws (`n` is an alias for `size`).

- obs:

  Observation to condition on (defaults to the posterior's `x_obs`).

- refresh:

  Force a new run even when a cached one would do.

- verbose:

  Report sampling progress.

- ...:

  Unused.

## Value

An `n x dim` matrix of posterior draws (class `nsbi_samples`), with
convergence diagnostics attached as attribute `diagnostics`.
