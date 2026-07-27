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

Numeric vector: the MAP parameter estimate.

## Details

On a posterior from
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) the
initial draws come from MCMC, so `n_init` buys a chain rather than a
forward pass. They are cached on the posterior like any other run.
