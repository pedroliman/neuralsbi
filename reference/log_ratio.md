# Evaluate a learned likelihood-to-evidence ratio

`log_ratio()` evaluates the ratio learned by
[`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md), \\\log
r(\theta, x) = \log p(x \mid \theta) - \log p(x)\\.

## Usage

``` r
log_ratio(fit, theta, x, ...)

# S3 method for class 'nsbi_nre'
log_ratio(fit, theta, x, sum_iid = TRUE, max_batch = 1e+05, ...)
```

## Arguments

- fit:

  An `nsbi_nre` fit from
  [`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md).

- theta:

  Parameter values: a numeric vector (one parameter set) or an
  `n_theta x dim_theta` matrix.

- x:

  Observed data: a numeric vector (one observation) or an
  `n_obs x dim_x` matrix whose rows are independent observations.

- ...:

  Unused, for S3 consistency.

- sum_iid:

  Sum the log ratio over the rows of `x` (the default). Set `FALSE` to
  get the per-observation values instead.

- max_batch:

  Largest number of `(theta, x)` pairs evaluated in one call to the
  classifier. Only affects memory and speed.

## Value

With `sum_iid = TRUE`, a numeric vector with one entry per row of
`theta`. With `sum_iid = FALSE`, an `n_theta x n_obs` matrix.

## Details

Rows of `x` are treated as **independent observations from the same
parameter**, so by default the result sums over them. That sum is the
log-likelihood of the whole data set up to an additive constant that
does not depend on `theta`, which is all a posterior needs and all NRE
learns: the evidence term \\\log p(x)\\ is never estimated separately,
so `log_ratio()` is not a log-likelihood you can compare across
observations. Differences between two `theta` at the same `x` are
exactly the differences in log-likelihood.

## See also

[`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md)
to turn the ratio into posterior draws,
[`log_lik()`](https://neuralsbi.pedrodelima.com/reference/log_lik.md)
for the [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md)
equivalent that *is* a calibrated log-likelihood.

## Examples

``` r
prior <- prior_uniform(c(mu = -3), c(mu = 3))
fit <- nre(prior, function(mu) c(y = rnorm(1, mu, 0.5)),
           n_simulations = 1000, classifier = "logistic")

x_obs <- matrix(rnorm(20, mean = 1, sd = 0.5), ncol = 1)
grid <- matrix(seq(-2, 2, length.out = 5), ncol = 1)
log_ratio(fit, grid, x_obs)
#> [1] -370.224203 -141.967041    3.034027   64.779000   43.267879
```
