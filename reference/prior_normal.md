# Independent normal prior

Independent normal prior

## Usage

``` r
prior_normal(mean, sd = 1)
```

## Arguments

- mean:

  Numeric vector of means (one per parameter).

- sd:

  Numeric scalar or vector of standard deviations.

## Value

An `nsbi_prior` object.

## Examples

``` r
prior <- prior_normal(mean = c(0, 0), sd = 1)
```
