# Posterior from a neural likelihood

Bayes' rule turns an
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) fit into a
posterior, \\p(\theta \mid x) \propto q\_\phi(x \mid
\theta)\\p(\theta)\\, but the result has no closed form and no direct
sampler.
[`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md)
on an `nsbi_nle` therefore returns an object that samples with MCMC, and
[`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md) on
it runs a chain rather than a forward pass.

## Usage

``` r
# S3 method for class 'nsbi_nle'
posterior(
  fit,
  x_obs = NULL,
  sampler = c("slice", "stan"),
  n_chains = NULL,
  warmup = 200L,
  thin = 2L,
  init_strategy = c("resample", "proposal"),
  seed = NULL,
  ...
)
```

## Arguments

- fit:

  An `nsbi_nle` object from
  [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md).

- x_obs:

  Observation to condition on. Rows are treated as independent
  observations from the same parameter, and the log-likelihood sums over
  them.

- sampler:

  `"slice"` (the default) or `"stan"`.

- n_chains:

  Number of chains. The default depends on the sampler: 20 for
  `"slice"`, which evaluates every chain in one batched call per step
  and so pays almost nothing for more of them, and 4 for `"stan"`, where
  each chain is a separate process with its own warmup to pay for.

- warmup:

  Steps discarded at the start of each chain.

- thin:

  Keep one draw in `thin`. The default is 2: the slice width is adapted
  during warmup, and with that the retained draws are already close to
  independent – on a Gaussian target with 20 chains, `thin = 2` gives a
  bulk ESS of about 96% of them. Since each evaluation is a forward pass
  over every observation, thinning harder buys very little for what it
  costs. Raise it if the reported ESS says you need to. (Python `sbi`
  thins by 1 by default; it thinned by 10 up to v0.21.)

- init_strategy:

  `"resample"` (default, and `sbi`'s) weights a pool of prior draws by
  the posterior density and resamples the starting points from it,
  drawing more pools as needed if the first one lands fewer than
  `n_chains` draws in the posterior's support; `"proposal"` skips the
  weighting and keeps whichever prior draws land inside the posterior's
  support, also drawing more pools as needed. `"proposal"` is cheaper
  per accepted draw, but every draw the posterior excludes is wasted: if
  only a fraction `a` of the prior's mass survives, roughly `1 / a`
  draws are needed per starting point, so a posterior that rules out
  most of the prior can make `"proposal"` far slower than `"resample"`
  to find enough starting points, and it errors out once its draw budget
  is spent. The pool is 1000 draws, set with `n_pool`. `sbi` uses 10,000
  *per chain*, which on a surrogate summed over thousands of
  observations is a large bill before the first MCMC step.

- seed:

  Optional integer seed.

- ...:

  Further arguments to the sampler: `width`, `max_steps` and `n_pool`
  for `"slice"`, or `iter_warmup`, `iter_sampling` and `refresh` for
  `"stan"`.

## Value

An object of class
`c("nsbi_nle_posterior", "nsbi_mcmc_posterior", "nsbi_posterior")`.

## Details

Two samplers are available. `"slice"` is the default: a vectorized
univariate slice sampler (see
[nsbi_mcmc](https://neuralsbi.pedrodelima.com/reference/nsbi_mcmc.md))
with nothing to tune and no dependencies. `"stan"` writes the fitted
likelihood out as a Stan program (see
[`stan_code()`](https://neuralsbi.pedrodelima.com/reference/stan_export.md))
and runs NUTS on it through cmdstanr or rstan, which mixes better on
correlated posteriors at the cost of a one-time model compile.

## Examples

``` r
prior <- prior_uniform(c(mu = -3), c(mu = 3))
fit <- nle(prior, function(mu) c(y = rnorm(1, mu, 0.5)),
           n_simulations = 1000, density_estimator = "linear_gaussian")

x_obs <- matrix(rnorm(50, mean = 1, sd = 0.5), ncol = 1)
post <- posterior(fit, x_obs, n_chains = 4, warmup = 50, thin = 2)
draws <- sample(post, 400)
colMeans(draws)
#>       mu 
#> 1.186398 
```
