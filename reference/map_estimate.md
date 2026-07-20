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
