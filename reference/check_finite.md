# Require every entry to be finite

`NA`, `NaN` and `Inf` all reach the estimators as a
[`chol()`](https://rdrr.io/r/base/chol.html) failure or a non-finite
validation loss, which blames training for a bad input. Naming the
argument and the first offending position at the boundary is cheaper to
act on.

## Usage

``` r
check_finite(m, arg)
```

## Arguments

- m:

  A numeric vector or matrix.

- arg:

  Name of the argument.

## Value

`m`, invisibly.
