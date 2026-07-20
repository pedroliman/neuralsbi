# Box-uniform (independent uniform) prior

Box-uniform (independent uniform) prior

## Usage

``` r
prior_uniform(low, high)
```

## Arguments

- low:

  Numeric vector of lower bounds (one per parameter).

- high:

  Numeric vector of upper bounds (one per parameter).

## Value

An `nsbi_prior` object.

## Examples

``` r
prior <- prior_uniform(low = c(-2, -2, -2), high = c(2, 2, 2))
theta <- sample_prior(prior, 5)
```
