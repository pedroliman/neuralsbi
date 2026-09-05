# Get (or build and cache) the surrogate potential behind an MCMC posterior

[`surrogate_potential()`](https://neuralsbi.pedrodelima.com/reference/surrogate_potential.md)
standardizes `x_obs` and lets the estimator precompute whatever it can
from the observation alone (see its own docs) – work meant to be paid
once per observation, not once per evaluation.
[`slice_sample_surrogate()`](https://neuralsbi.pedrodelima.com/reference/slice_sample_surrogate.md)
already gets this for free because it builds the closure once and calls
it for every MCMC step of a whole chain, but
[`mcmc_log_prob()`](https://neuralsbi.pedrodelima.com/reference/mcmc_log_prob.md)
used to rebuild the closure on every single
[`log_prob()`](https://neuralsbi.pedrodelima.com/reference/log_prob.md)
call. That is expensive on its own, and ruinous for
[`map_estimate()`](https://neuralsbi.pedrodelima.com/reference/map_estimate.md),
whose optimizer calls
[`log_prob()`](https://neuralsbi.pedrodelima.com/reference/log_prob.md)
once per evaluation against the same `x_obs` (#272). The closure is
cached here on the posterior's cache environment and only rebuilt when
`x_obs` or `max_batch` changes.

## Usage

``` r
cached_surrogate_potential(post, x_obs, max_batch)
```

## Arguments

- post:

  An `nsbi_mcmc_posterior` object.

- x_obs:

  The (already resolved) observation to condition on.

- max_batch:

  Largest batch size to pass to
  [`surrogate_potential()`](https://neuralsbi.pedrodelima.com/reference/surrogate_potential.md).

## Value

`function(theta)`, the cached (or freshly built) potential.
