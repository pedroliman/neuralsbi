# Log ratio of a fitted ratio estimator

The ratio estimators' counterpart to
[de_log_prob()](https://neuralsbi.pedrodelima.com/reference/density_estimator.md):
one number per row, \\\log r(\theta, x)\\, in standardized space (which
for a ratio is also the original space – see
[`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md)). A
one-row `x` is broadcast against a taller `theta`, as `de_log_prob()`
does.

## Usage

``` r
de_log_ratio(de, theta, x)
```

## Arguments

- de:

  A fitted ratio estimator.

- theta:

  Standardized parameters.

- x:

  Standardized data.

## Value

A numeric vector with one entry per row.
