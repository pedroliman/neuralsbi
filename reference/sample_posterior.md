# Sample from a posterior (non-generic alias)

Identical to `sample(post, n)`; provided for users who prefer not to
rely on the generic.

## Usage

``` r
sample_posterior(post, n = 1000, obs = NULL, ...)
```

## Arguments

- post:

  An `nsbi_posterior` object.

- n:

  Number of posterior draws.

- obs:

  Observation to condition on (defaults to the posterior's `x_obs`).

- ...:

  Passed to
  [`sample.nsbi_posterior()`](https://pedroliman.github.io/neuralsbi/reference/sample.nsbi_posterior.md).

## Value

An `n x dim` matrix of posterior draws.
