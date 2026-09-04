# A surrogate likelihood as a plain R function

Fixes an observation and returns `function(theta)`, so the learned
likelihood can be handed to anything in R that wants a log-density:
[`stats::optim()`](https://rdrr.io/r/stats/optim.html), an MCMC package,
an importance sampler, a profile likelihood. The returned function is
vectorized – pass a matrix of parameter values and get one
log-likelihood per row.

## Usage

``` r
likelihood_fn(fit, x_obs, ...)
```

## Arguments

- fit:

  An `nsbi_nle` fit from
  [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md).

- x_obs:

  The observation to condition on. Rows are independent observations;
  the log-likelihood sums over them.

- ...:

  Passed to
  [`log_lik()`](https://neuralsbi.pedrodelima.com/reference/log_lik.md).

## Value

A function of `theta` returning a numeric vector of log-likelihoods,
with the observation attached as attribute `x_obs`.

## Examples

``` r
prior <- prior_uniform(c(mu = -3), c(mu = 3))
fit <- nle(prior, function(mu) c(y = rnorm(1, mu, 0.5)),
           n_simulations = 1000, density_estimator = "linear_gaussian")
#> Running the simulator sequentially. To use all your cores:
#>   library(future)
#>   plan(multisession)
#> Hide this hint with options(neuralsbi.parallel_hint = FALSE).

loglik <- likelihood_fn(fit, matrix(rnorm(20, 1, 0.5), ncol = 1))
optimize(function(m) loglik(m), c(-3, 3), maximum = TRUE)$maximum
#> [1] 0.9550394
```
