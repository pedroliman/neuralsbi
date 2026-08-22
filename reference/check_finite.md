# Require every entry to be finite

`NA`, `NaN` and `Inf` all reach the estimators as a
[`chol()`](https://rdrr.io/r/base/chol.html) failure or a non-finite
validation loss, which blames training for a bad input. Naming the
argument and the first offending position at the boundary is cheaper to
act on.

## Usage

``` r
check_finite(m, arg, allow_inf = FALSE)
```

## Arguments

- m:

  A numeric vector or matrix.

- arg:

  Name of the argument.

- allow_inf:

  Let `Inf`/`-Inf` through and only reject `NA`/`NaN`. Use this where
  `Inf` already has a well-defined meaning downstream: an MCMC
  posterior's `theta` routes through the prior's own density first,
  which correctly sends an infinite parameter value to zero mass rather
  than to the estimator.

## Value

`m`, invisibly.
