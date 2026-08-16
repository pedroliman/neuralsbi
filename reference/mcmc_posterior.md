# Check the arguments of an MCMC-sampled posterior and assemble it

Shared by
[`posterior.nsbi_nle()`](https://neuralsbi.pedrodelima.com/reference/posterior.nsbi_nle.md)
and
[`posterior.nsbi_nre()`](https://neuralsbi.pedrodelima.com/reference/posterior.nsbi_nre.md).
Both wrap a fit and a sampler configuration around a draw cache, and
every argument except the sampler is checked the same way; only the
class they carry and the samplers they allow differ, and both of those
are settled by the caller before it gets here.

## Usage

``` r
mcmc_posterior(
  fit,
  x_obs,
  sampler,
  n_chains,
  warmup,
  thin,
  init_strategy,
  seed,
  dots,
  class
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
  the posterior density and resamples the starting points from it;
  `"proposal"` skips the weighting and keeps whichever prior draws land
  inside the posterior's support, drawing more pools as needed.
  `"proposal"` is cheaper per accepted draw, but every draw the
  posterior excludes is wasted: if only a fraction `a` of the prior's
  mass survives, roughly `1 / a` draws are needed per starting point, so
  a posterior that rules out most of the prior can make `"proposal"` far
  slower than `"resample"` to find enough starting points, and it errors
  out once its draw budget is spent. The pool is 1000 draws, set with
  `n_pool`. `sbi` uses 10,000 *per chain*, which on a surrogate summed
  over thousands of observations is a large bill before the first MCMC
  step.

- seed:

  Optional integer seed.

- dots:

  The sampler arguments the caller collected from `...`.

- class:

  The posterior class to stamp on the result.
