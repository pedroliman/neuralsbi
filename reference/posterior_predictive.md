# Posterior predictive draws

Samples parameters from the posterior and pushes them back through the
simulator, giving predictive data to compare against the observation. A
draw whose simulation returns non-finite output is dropped, which lowers
the number of predictive draws returned.

## Usage

``` r
posterior_predictive(post, simulator, n = 1000L, x = NULL, sim_args = list())
```

## Arguments

- post:

  An `nsbi_posterior` object.

- simulator:

  The simulator; called once per posterior draw (see
  [nsbi_simulator](https://neuralsbi.pedrodelima.com/reference/nsbi_simulator.md)).

- n:

  Number of predictive draws.

- x:

  Observation to condition on (defaults to `x_obs`).

- sim_args:

  Named list of extra arguments passed to every simulator call; see
  [nsbi_simulator](https://neuralsbi.pedrodelima.com/reference/nsbi_simulator.md).

## Value

An `n x d` matrix of simulated data from posterior parameter draws.
