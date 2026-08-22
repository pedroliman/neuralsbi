# Evaluate a surrogate likelihood

`log_lik()` evaluates the likelihood learned by
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md): \\\log
q\_\phi(x \mid \theta)\\, in the original units of both.

## Usage

``` r
log_lik(fit, theta, x, ...)

# S3 method for class 'nsbi_nle'
log_lik(fit, theta, x, sum_iid = TRUE, max_batch = 1e+05, ...)
```

## Arguments

- fit:

  An `nsbi_nle` fit from
  [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md).

- theta:

  Parameter values: a numeric vector (one parameter set) or an
  `n_theta x dim_theta` matrix.

- x:

  Observed data: a numeric vector (one observation) or an
  `n_obs x dim_x` matrix whose rows are independent observations.

- ...:

  Unused, for S3 consistency.

- sum_iid:

  Sum the log-density over the rows of `x` (the default). Set `FALSE` to
  get the per-observation values instead.

- max_batch:

  Largest number of `(theta, x)` pairs evaluated in one call to the
  estimator. Only affects memory and speed.

## Value

With `sum_iid = TRUE`, a numeric vector with one entry per row of
`theta`. With `sum_iid = FALSE`, an `n_theta x n_obs` matrix.

## Details

Rows of `x` are treated as **independent observations from the same
parameter**, so by default the result sums over them, \\\sum_i \log
q\_\phi(x_i \mid \theta)\\. That sum is the whole point of NLE: the
estimator is trained on one observation at a time, and the number of
observations you condition on afterwards is free.

## See also

[`likelihood_fn()`](https://neuralsbi.pedrodelima.com/reference/likelihood_fn.md)
for a closure over a fixed observation,
[`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md)
to turn the likelihood into posterior draws.

## Examples

``` r
prior <- prior_uniform(c(mu = -3), c(mu = 3))
fit <- nle(prior, function(mu) c(y = rnorm(1, mu, 0.5)),
           n_simulations = 1000, density_estimator = "linear_gaussian")

x_obs <- matrix(rnorm(20, mean = 1, sd = 0.5), ncol = 1)
grid <- matrix(seq(-2, 2, length.out = 5), ncol = 1)
log_lik(fit, grid, x_obs)
#> [1] -363.40767 -168.41256  -52.41733  -15.42197  -57.42650
```
