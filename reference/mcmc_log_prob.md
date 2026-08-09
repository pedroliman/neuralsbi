# The unnormalized log density behind an MCMC posterior's [`log_prob()`](https://neuralsbi.pedrodelima.com/reference/log_prob.md) method

Neither surrogate gives the evidence \\p(x)\\, so `normalize = TRUE` has
nothing to normalize by and says so rather than returning a number that
looks like a density.

## Usage

``` r
mcmc_log_prob(post, theta, x, warn, what)
```

## Arguments

- post:

  An MCMC-sampled `nsbi_posterior`.

- theta:

  Parameter values to evaluate.

- x:

  Observation to condition on, or `NULL` for the posterior's own.

- warn:

  Warn that `normalize` is being ignored (the caller decides, since only
  it can see whether the argument was actually supplied).

- what:

  The method name to use in that warning, `"NLE"` or `"NRE"`.
