# Box-uniform (independent uniform) prior

Box-uniform (independent uniform) prior

## Usage

``` r
prior_uniform(low, high)
```

## Arguments

- low:

  Numeric vector of lower bounds (one per parameter). Naming the vector
  (e.g. `c(beta = 0, gamma = 0)`) attaches those names to every
  downstream parameter matrix, posterior sample, and diagnostic plot.

- high:

  Numeric vector of upper bounds (one per parameter).

## Value

An `nsbi_prior` object.

## Examples

``` r
prior <- prior_uniform(low = c(-2, -2, -2), high = c(2, 2, 2))
theta <- sample_prior(prior, 5)
```
