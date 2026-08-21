# Location-scale Student-t

[`stats::dt()`](https://rdrr.io/r/stats/TDist.html) only knows the
standardized t, and Stan's `student_t` is the location-scale one. These
three wrap the shift and the scale, including the `-log(scale)` Jacobian
on the density, and take their arguments under the names the registry
uses.

## Usage

``` r
dstudent_t(x, df, location = 0, scale = 1, log = FALSE)

pstudent_t(q, df, location = 0, scale = 1)

qstudent_t(p, df, location = 0, scale = 1)
```

## Arguments

- x, q, p:

  Value, quantile or probability.

- df, location, scale:

  Distribution parameters.

- log:

  Return the log density.
