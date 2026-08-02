# Warn about columns standardization cannot scale

Names the columns by name where they have one and by index otherwise,
since an index is no help when the matrix came from a data frame and a
name is not available when it did not. This is a warning rather than an
error because a single row is a legitimate way to get here:
[`sd()`](https://rdrr.io/r/stats/sd.html) of one value is `NA`.

## Usage

``` r
warn_constant_columns(x, flat, what)
```

## Arguments

- x:

  The matrix being standardized.

- flat:

  Logical vector marking the columns with no usable spread.

- what:

  Name of the argument `x` came from, `"theta"` or `"x"`.
