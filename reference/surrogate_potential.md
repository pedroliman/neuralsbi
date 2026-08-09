# Unnormalized log posterior of a surrogate fit

\\\log q\_\phi(x \mid \theta) + \log p(\theta)\\ for an
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) fit and
\\\log r\_\phi(\theta, x) + \log p(\theta)\\ for an
[`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md) fit,
returning `-Inf` outside the prior support. This is the potential the
MCMC samplers target, and the two fits differ only in what
[`surrogate_ops()`](https://neuralsbi.pedrodelima.com/reference/surrogate_ops.md)
hands back.

## Usage

``` r
surrogate_potential(fit, x_obs, max_batch = 1e+05)
```

## Arguments

- fit:

  An `nsbi_nle` or `nsbi_nre` fit.

- x_obs:

  The observation to condition on. Rows are independent observations of
  the same parameter.

- max_batch:

  Largest number of `(theta, x)` pairs evaluated at once. No caller
  overrides it today, but it is a real tuning knob – the batch size the
  MCMC evaluations are chunked into – that a future caller would
  plausibly want to change, so it stays a parameter.

## Value

`function(theta)` giving one unnormalized log posterior density per row.
