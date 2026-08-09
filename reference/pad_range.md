# Data range, padded like a default ggplot2 continuous scale

[`pairplot()`](https://neuralsbi.pedrodelima.com/reference/pairplot.md)'s
default axis limits (when the caller does not pass `limits`) come from
this: the raw [`range()`](https://rdrr.io/r/base/range.html) of a
column, expanded 5% on each side, matching ggplot2's own default
continuous-scale expansion so the shared range still leaves a little air
around the extreme draws. A constant column (`diff(r) == 0`) would
otherwise pad by nothing, so that case gets a small absolute pad
instead.

## Usage

``` r
pad_range(x)
```

## Arguments

- x:

  Numeric vector.

## Value

`c(lo, hi)`.
