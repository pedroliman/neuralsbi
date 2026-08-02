# Learn a standardizer from a matrix

A column with no spread cannot be divided by its standard deviation, so
it keeps scale 1. Pass `what` to hear about it: the guard is silent
otherwise, and a constant column is worth a word because nothing
downstream will complain. Training converges, the posterior looks
plausible, and the coordinate does nothing.

## Usage

``` r
fit_standardizer(x, eps = 1e-08, what = NULL)
```

## Arguments

- x:

  The matrix to learn from.

- eps:

  Standard deviations below this count as no spread.

- what:

  Name of the argument `x` came from (`"theta"` or `"x"`), used in the
  warning. `NULL`, the default, warns about nothing. The
  `standardize = FALSE` path in
  [`prepare_simulations()`](https://neuralsbi.pedrodelima.com/reference/prepare_simulations.md)
  builds a degenerate standardizer from a one-row zero matrix on
  purpose, and the diagnostics standardize draws they generated
  themselves.
