# Posterior objects

A posterior wraps a trained
[`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) fit
together with (optionally) a default observation `x_obs`. It knows how
to draw posterior samples, evaluate the posterior log-density, and find
the maximum-a-posteriori (MAP) estimate. All transforms between
standardized training space and the original parameter space are handled
internally.

The three inference methods reach a posterior by different routes, and
`posterior()` hides the difference. An
[`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) fit
already *is* a posterior estimator, so the returned object samples with
a forward pass. An
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) fit only
knows the likelihood and an
[`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md) fit only
the likelihood ratio, so both return an object that samples with MCMC
and takes the extra arguments that implies.

## Usage

``` r
posterior(fit, x_obs = NULL, ...)

# Default S3 method
posterior(fit, x_obs = NULL, ...)

# S3 method for class 'nsbi_npe'
posterior(fit, x_obs = NULL, ...)
```

## Arguments

- fit:

  An `nsbi_npe` object from
  [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) or
  [`npe_sequential()`](https://neuralsbi.pedrodelima.com/reference/npe_sequential.md),
  an `nsbi_nle` object from
  [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md), or an
  `nsbi_nre` object from
  [`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md).

- x_obs:

  Optional default observation to condition on. If supplied it becomes
  the default `x` for
  [`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md),
  [`log_prob()`](https://neuralsbi.pedrodelima.com/reference/log_prob.md)
  and
  [`map_estimate()`](https://neuralsbi.pedrodelima.com/reference/map_estimate.md).
  For an NLE or NRE fit, rows of `x_obs` are independent observations.

- ...:

  Passed to methods. See
  [`posterior.nsbi_nle()`](https://neuralsbi.pedrodelima.com/reference/posterior.nsbi_nle.md)
  and
  [`posterior.nsbi_nre()`](https://neuralsbi.pedrodelima.com/reference/posterior.nsbi_nre.md)
  for the MCMC controls those fits accept.

## Value

An `nsbi_posterior` object.

## Details

For bounded priors, samples that fall outside the prior support are
rejected ("leakage" correction), and
[`log_prob()`](https://neuralsbi.pedrodelima.com/reference/log_prob.md)
is renormalized by the estimated acceptance probability so it integrates
to one over the support.

## See also

[`save_npe()`](https://neuralsbi.pedrodelima.com/reference/save_npe.md),
which is how a torch-backed fit gets to disk and back;
[`readRDS()`](https://rdrr.io/r/base/readRDS.html) returns one whose
network is dead, and `posterior()` says so rather than failing later.
