# Named prior families

Priors built from a named distribution family, one marginal per
parameter, following Stan's argument names and argument order. Each
constructor is vectorized: pass a vector where Stan would write one
sampling statement per component, and pass a scalar to use the same
value for every parameter. Naming any of the vectors (e.g.
`c(beta = 2, gamma = 3)`) names the parameters, exactly as naming `low`
does in
[`prior_uniform()`](https://neuralsbi.pedrodelima.com/reference/prior_uniform.md).

## Usage

``` r
prior_lognormal(meanlog = 0, sdlog = 1)

prior_exponential(rate = 1)

prior_gamma(shape, rate = 1)

prior_beta(shape1, shape2)

prior_student_t(df, location = 0, scale = 1)

prior_cauchy(location = 0, scale = 1)

prior_half_normal(sd = 1)

prior_half_cauchy(scale = 1)
```

## Arguments

- meanlog, sdlog:

  Log-scale mean and standard deviation of the log-normal, as in
  [`stats::dlnorm()`](https://rdrr.io/r/stats/Lognormal.html) and Stan's
  `lognormal`.

- rate:

  Rate of the exponential or gamma, as in
  [`stats::dexp()`](https://rdrr.io/r/stats/Exponential.html) and Stan's
  `exponential`/`gamma`. Must be positive.

- shape:

  Shape of the gamma. Must be positive.

- shape1, shape2:

  Beta shape parameters, as in
  [`stats::dbeta()`](https://rdrr.io/r/stats/Beta.html). Both must be
  positive.

- df:

  Degrees of freedom of the Student-t. Must be positive.

- location, scale:

  Location and scale of the Student-t, the Cauchy or the half-Cauchy.
  `scale` must be positive.

- sd:

  Standard deviation of the normal underlying the half-normal. Must be
  positive.

## Value

An `nsbi_prior` object.

## Details

Families with constrained support set `lower`/`upper` on the prior, so
the posterior's leakage correction (see
[`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md))
rejects and renormalizes against the right region without any further
declaration. A gamma prior bounds at zero, a beta prior at zero and one,
and the half families at zero.

A prior built here can be truncated further with
[`prior_truncated()`](https://neuralsbi.pedrodelima.com/reference/prior_truncated.md),
combined with others through
[`prior_independent()`](https://neuralsbi.pedrodelima.com/reference/prior_independent.md),
and written out by
[`stan_code()`](https://neuralsbi.pedrodelima.com/reference/stan_export.md)
as the sampling statement it came from. That is the reason these exist
as families rather than as recipes for
[`prior_custom()`](https://neuralsbi.pedrodelima.com/reference/prior_custom.md):
the family is what carries over to Stan and to the truncation constant.

`prior_half_normal()` and `prior_half_cauchy()` are the zero-centred
families truncated to \\\[0, \infty)\\, matching Stan's idiom of
declaring `real<lower=0>` and writing `sigma ~ normal(0, 1)`. Their
densities carry the `log(2)` renormalization, which Stan drops as a
constant and this package cannot. For a half-normal centred somewhere
other than zero, wrap
[`prior_normal()`](https://neuralsbi.pedrodelima.com/reference/prior_normal.md)
in
[`prior_truncated()`](https://neuralsbi.pedrodelima.com/reference/prior_truncated.md).

## See also

[`prior_independent()`](https://neuralsbi.pedrodelima.com/reference/prior_independent.md)
to combine several of these into one joint prior,
[`prior_truncated()`](https://neuralsbi.pedrodelima.com/reference/prior_truncated.md)
to bound one, and
[priors](https://neuralsbi.pedrodelima.com/reference/priors.md) for the
rest.

## Examples

``` r
# One log-normal per parameter, named through the vector.
prior <- prior_lognormal(meanlog = c(beta = log(0.4), gamma = log(0.125)),
                         sdlog = c(0.5, 0.2))
theta <- sample_prior(prior, 5)

# A scalar argument is shared: three gammas with the same rate.
prior_gamma(shape = c(2, 5, 9), rate = 3)
#> <nsbi_prior> type=gamma, dim=3
#>   lower: 0, 0, 0 
#>   theta[1] ~ gamma(2, 3)
#>   theta[2] ~ gamma(5, 3)
#>   theta[3] ~ gamma(9, 3)

# Support bounds come from the family, and drive the leakage correction.
within_support(prior_beta(2, 15), c(0.1, 1.5))
#> [1]  TRUE FALSE
```
