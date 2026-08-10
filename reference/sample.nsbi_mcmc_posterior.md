# Sample an MCMC posterior

Runs the sampler chosen by
[`posterior.nsbi_nle()`](https://neuralsbi.pedrodelima.com/reference/posterior.nsbi_nle.md)
or
[`posterior.nsbi_nre()`](https://neuralsbi.pedrodelima.com/reference/posterior.nsbi_nre.md).
Draws are cached on the posterior object, so asking for the same or
fewer draws from the same observation returns immediately instead of
re-running a chain – which is what makes
[`summary()`](https://rdrr.io/r/base/summary.html) and repeated calls
tolerable.

## Usage

``` r
# S3 method for class 'nsbi_mcmc_posterior'
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

  An `nsbi_mcmc_posterior` object (named `x` to satisfy the
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
