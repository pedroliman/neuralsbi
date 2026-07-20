# Posterior objects

A posterior wraps a trained
[`npe()`](https://pedroliman.github.io/neural.sbi/reference/npe.md) fit
together with (optionally) a default observation `x_obs`. It knows how
to draw posterior samples, evaluate the posterior log-density, and find
the maximum-a-posteriori (MAP) estimate. All transforms between
standardized training space and the original parameter space are handled
internally.

## Usage

``` r
posterior(fit, x_obs = NULL)
```

## Arguments

- fit:

  An `nsbi_npe` object from
  [`npe()`](https://pedroliman.github.io/neural.sbi/reference/npe.md).

- x_obs:

  Optional default observation to condition on. If supplied it becomes
  the default `x` for
  [`sample()`](https://pedroliman.github.io/neural.sbi/reference/sample.md),
  [`log_prob()`](https://pedroliman.github.io/neural.sbi/reference/log_prob.md)
  and
  [`map_estimate()`](https://pedroliman.github.io/neural.sbi/reference/map_estimate.md).

## Value

An `nsbi_posterior` object.

## Details

For bounded priors, samples that fall outside the prior support are
rejected ("leakage" correction), and
[`log_prob()`](https://pedroliman.github.io/neural.sbi/reference/log_prob.md)
is renormalized by the estimated acceptance probability so it integrates
to one over the support.
