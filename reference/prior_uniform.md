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
  `NA`/`NaN` are rejected; `-Inf` is accepted (matching `high`'s `Inf`)
  but makes the prior improper, so
  [`sample_prior()`](https://neuralsbi.pedrodelima.com/reference/sample_prior.md)
  on it errors unless it is first bounded with
  [`prior_truncated()`](https://neuralsbi.pedrodelima.com/reference/prior_truncated.md).

- high:

  Numeric vector of upper bounds (one per parameter). Same finiteness
  rule as `low`, with `Inf` in place of `-Inf`.

## Value

An `nsbi_prior` object.

## Examples

``` r
prior <- prior_uniform(low = c(-2, -2, -2), high = c(2, 2, 2))
theta <- sample_prior(prior, 5)
```
