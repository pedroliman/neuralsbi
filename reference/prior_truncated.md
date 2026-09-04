# Truncate a prior to a box

Stan's `T[lower, upper]`, as an object. The returned prior is the
original restricted to the box and renormalized by the mass it keeps, so
its `log_prob` is a proper log density rather than the original shifted
by an unknown constant. That matters here in a way it does not in Stan:
the density is compared against a learned posterior in
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md)'s MCMC
target and summed with it in
[`c2st()`](https://neuralsbi.pedrodelima.com/reference/c2st.md) and the
diagnostics, so a missing constant is a wrong answer rather than a
constant offset.

## Usage

``` r
prior_truncated(prior, lower = NULL, upper = NULL)
```

## Arguments

- prior:

  An `nsbi_prior` from a named family.

- lower, upper:

  Truncation bounds. Numeric of length `prior$dim`, or length 1 to apply
  the same bound to every parameter. Give one or both; `-Inf`/`Inf`
  leaves that side alone.

## Value

An `nsbi_prior` object.

## Details

Renormalization needs the family's CDF, so `prior` has to come from a
named family:
[`prior_uniform()`](https://neuralsbi.pedrodelima.com/reference/prior_uniform.md),
[`prior_normal()`](https://neuralsbi.pedrodelima.com/reference/prior_normal.md),
anything in
[prior_families](https://neuralsbi.pedrodelima.com/reference/prior_families.md),
or a
[`prior_independent()`](https://neuralsbi.pedrodelima.com/reference/prior_independent.md)
built from those. A
[`prior_custom()`](https://neuralsbi.pedrodelima.com/reference/prior_custom.md)
is arbitrary R code with no CDF behind it; give it `lower`/`upper` there
instead, which rejects out-of-support draws without claiming to
renormalize.

## See also

[prior_families](https://neuralsbi.pedrodelima.com/reference/prior_families.md),
[`prior_independent()`](https://neuralsbi.pedrodelima.com/reference/prior_independent.md),
[priors](https://neuralsbi.pedrodelima.com/reference/priors.md).

## Examples

``` r
# A half-normal, the long way round (prior_half_normal() is the short one).
prior_truncated(prior_normal(mean = 0, sd = 1), lower = 0)
#> <nsbi_prior> type=truncated, dim=1
#>   lower: 0 
#>   theta[1] ~ normal(0, 1) T[0, ]

# A contact rate known to be between 0.1 and 2, log-normal in between.
prior <- prior_truncated(prior_lognormal(log(0.4), 0.5),
                         lower = 0.1, upper = 2)
range(sample_prior(prior, 100))
#> [1] 0.1397014 1.4725125
```
