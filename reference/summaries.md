# Summaries and tidy accessors

[`summary()`](https://rdrr.io/r/base/summary.html) methods for fits,
posteriors, and samples, plus
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) for
posterior draws so results drop straight into data-frame workflows
(dplyr, ggplot2, ...).

## Usage

``` r
# S3 method for class 'nsbi_samples'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'nsbi_samples'
summary(object, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), ...)

# S3 method for class 'nsbi_posterior'
summary(object, n = 1000L, x = NULL, ...)

# S3 method for class 'nsbi_npe'
summary(object, ...)
```

## Arguments

- x:

  Observation to condition on (defaults to the posterior's `x_obs`); for
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), the
  samples object.

- row.names, optional:

  Standard
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  arguments.

- ...:

  Additional arguments passed to methods.

- object:

  An `nsbi_samples` matrix from
  [`sample()`](https://pedroliman.github.io/neuralsbi/reference/sample.md),
  an `nsbi_posterior`, or an `nsbi_npe` fit.

- probs:

  Quantiles to report.

- n:

  Number of draws used to summarize a posterior.

## Value

For samples and posteriors, a data frame with one row per parameter
(mean, sd, and quantiles). For fits, an invisible list of training
metadata.
