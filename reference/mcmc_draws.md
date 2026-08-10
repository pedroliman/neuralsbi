# Run (or reuse) the chain behind an MCMC posterior's [`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md) method

The body
[`sample.nsbi_mcmc_posterior()`](https://neuralsbi.pedrodelima.com/reference/sample.nsbi_mcmc_posterior.md)
runs regardless of whether the underlying fit is an
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) or
[`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md). Which
sampler runs is read off the posterior object, so nothing here needs to
know which kind of fit produced it.

## Usage

``` r
mcmc_draws(post, n, obs, refresh, verbose)
```

## Arguments

- post:

  An MCMC-sampled `nsbi_posterior`.

- n:

  Number of draws.

- obs:

  Observation to condition on, or `NULL` for the posterior's own.

- refresh:

  Force a new run even when a cached one would do.

- verbose:

  Report sampling progress.
