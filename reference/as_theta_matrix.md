# Coerce parameters/data to a numeric matrix with a known column count

Preserves column names where they carry meaning: a data frame's or
matrix's existing `colnames`, or a plain vector's
[`names()`](https://rdrr.io/r/base/names.html) when the vector is
interpreted as a single row. A vector re-interpreted as a stacked column
(the `byrow` branch) loses its names – they described entries, not a
shared parameter/outcome identity.

## Usage

``` r
as_theta_matrix(x, d = NULL)
```
