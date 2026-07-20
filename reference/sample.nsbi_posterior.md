# Sample from a posterior

Sample from a posterior

## Usage

``` r
# S3 method for class 'nsbi_posterior'
sample(x, size = 1000, n = size, obs = NULL, max_sampling_batches = 100L, ...)
```

## Arguments

- x:

  An `nsbi_posterior` object (named `x` to satisfy the
  [`sample()`](https://pedroliman.github.io/neural.sbi/reference/sample.md)
  generic).

- size, n:

  Number of posterior draws (`n` is an alias for `size`).

- obs:

  Observation to condition on (defaults to the posterior's `x_obs`).

- max_sampling_batches:

  Safety cap on rejection-sampling rounds for bounded priors.

- ...:

  Unused.

## Value

An `n x dim` matrix of posterior draws (class `nsbi_samples`).
