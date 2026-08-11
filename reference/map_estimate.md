# Maximum a posteriori (MAP) estimate

Starts from the best of a set of posterior draws and refines with a
derivative-free optimizer.

## Usage

``` r
map_estimate(post, x = NULL, n_init = 1000L)
```

## Arguments

- post:

  An `nsbi_posterior` object.

- x:

  Observation to condition on (defaults to `x_obs`).

- n_init:

  Number of initial draws used to seed the search.

## Value

Numeric vector: the MAP parameter estimate. For a bounded prior (from
[`prior_uniform()`](https://neuralsbi.pedrodelima.com/reference/prior_uniform.md)
or a
[`prior_custom()`](https://neuralsbi.pedrodelima.com/reference/prior_custom.md)
with `lower`/`upper`), the estimate always falls inside the prior's
support – the search never accepts a step that leaves it, the same
guarantee
[`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md) and
[`log_prob()`](https://neuralsbi.pedrodelima.com/reference/log_prob.md)
give. Errors if the seeding draw comes back short of `n_init` –
including empty – which for a bounded prior means the estimator is
leaking mass outside the prior support faster than rejection sampling
can keep up; there is no starting point to search from in that case, so
this stops rather than continuing on a shorter, silently misleading
draw.

## Details

On a posterior from
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) or
[`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md) the
initial draws come from MCMC, so `n_init` buys a chain rather than a
forward pass. They are cached on the posterior like any other run.
