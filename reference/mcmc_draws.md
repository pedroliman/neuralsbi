# Run (or reuse) the chain behind an MCMC posterior's [`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md) method

The body
[`sample.nsbi_nle_posterior()`](https://neuralsbi.pedrodelima.com/reference/sample.nsbi_nle_posterior.md)
and
[`sample.nsbi_nre_posterior()`](https://neuralsbi.pedrodelima.com/reference/sample.nsbi_nre_posterior.md)
share. Which sampler runs is read off the posterior object, so nothing
here needs to know which kind of fit produced it.

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
