# Posterior predictive draws

Samples parameters from the posterior and pushes them back through the
simulator, giving predictive data to compare against the observation.

## Usage

``` r
posterior_predictive(post, simulator, n = 1000L, x = NULL)
```

## Arguments

- post:

  An `nsbi_posterior` object.

- simulator:

  The simulator.

- n:

  Number of predictive draws.

- x:

  Observation to condition on (defaults to `x_obs`).

## Value

An `n x d` matrix of simulated data from posterior parameter draws.
