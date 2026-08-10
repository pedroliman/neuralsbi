# Posterior from a neural ratio

Bayes' rule turns an
[`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md) fit into a
posterior, \\p(\theta \mid x) \propto r\_\phi(\theta, x)\\p(\theta)\\,
and as with
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) the result
has no closed form and no direct sampler.
[`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md)
on an `nsbi_nre` therefore returns an object that samples with MCMC, and
[`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md) on
it runs a chain rather than a forward pass.

## Usage

``` r
# S3 method for class 'nsbi_nre'
posterior(
  fit,
  x_obs = NULL,
  n_chains = 20L,
  warmup = 200L,
  thin = 2L,
  init_strategy = c("resample", "proposal"),
  seed = NULL,
  ...
)
```

## Arguments

- fit:

  An `nsbi_nre` object from
  [`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md).

- x_obs:

  Observation to condition on. Rows are treated as independent
  observations from the same parameter, and the log ratio sums over
  them.

- n_chains:

  Number of chains. The slice sampler evaluates every chain in one
  batched call per step, so more chains cost almost nothing.

- warmup:

  Steps discarded at the start of each chain.

- thin:

  Keep one draw in `thin`. See
  [`posterior.nsbi_nle()`](https://neuralsbi.pedrodelima.com/reference/posterior.nsbi_nle.md)
  for why the default is 2.

- init_strategy:

  `"resample"` (default, and `sbi`'s) weights a pool of prior draws by
  the posterior density and resamples the starting points from it;
  `"proposal"` starts from plain prior draws.

- seed:

  Optional integer seed.

- ...:

  Further arguments to the sampler: `width`, `max_steps` and `n_pool`.

## Value

An object of class
`c("nsbi_nre_posterior", "nsbi_mcmc_posterior", "nsbi_posterior")`.

## Details

There is one sampler here, the vectorized slice sampler of
[nsbi_mcmc](https://neuralsbi.pedrodelima.com/reference/nsbi_mcmc.md),
and no `sampler` argument to choose it with.
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md)'s `"stan"`
option works by transpiling the fitted density into a Stan `functions`
block; the residual classifier behind a ratio estimator has no such
export, so NUTS is not on offer.

## Examples

``` r
prior <- prior_uniform(c(mu = -3), c(mu = 3))
fit <- nre(prior, function(mu) c(y = rnorm(1, mu, 0.5)),
           n_simulations = 1000, classifier = "logistic")

x_obs <- matrix(rnorm(50, mean = 1, sd = 0.5), ncol = 1)
post <- posterior(fit, x_obs, n_chains = 4, warmup = 50, thin = 2)
draws <- sample(post, 400)
colMeans(draws)
#>        mu 
#> 0.9817372 
```
