# Plot posterior predictive checks

Compares data simulated from posterior parameter draws (see
[`posterior_predictive()`](https://pedroliman.github.io/neuralsbi/reference/posterior_predictive.md))
with the observed data, one marginal histogram per data dimension with
the observation marked. If the observation falls in the tails of the
predictive distribution, the model (or the fit) does not reproduce the
data it is conditioned on.

## Usage

``` r
plot_posterior_predictive(pred, x_obs, labels = NULL, bins = 30L)
```

## Arguments

- pred:

  A matrix of predictive draws from
  [`posterior_predictive()`](https://pedroliman.github.io/neuralsbi/reference/posterior_predictive.md).

- x_obs:

  The observed data vector the posterior was conditioned on.

- labels:

  Optional labels for the data dimensions.

- bins:

  Number of histogram bins.

## Value

A `ggplot` object (also drawn as a side effect), invisibly.
