# Conditional mean of the linear-Gaussian estimator

`lingauss_mean()` is the only place `x` meets `de$B`, so it is also the
place the width of `x` is checked. It passes `de$dim_x` to
[`as_theta_matrix()`](https://neuralsbi.pedrodelima.com/reference/as_theta_matrix.md)
for the reason every neural estimator passes its own: a wrong-width `x`
otherwise gets as far as the matrix product and is reported as
"non-conformable arguments", which names neither the argument nor the
width expected of it. An estimator fitted before `dim_x` was recorded
has `NULL` here and keeps the old unchecked behaviour.

## Usage

``` r
lingauss_mean(de, x)
```

## Arguments

- de:

  A fitted `nsbi_de_lingauss` object.

- x:

  The conditioning variable: a numeric vector, matrix or data frame with
  `de$dim_x` columns.

## Value

An `nrow(x) x de$dim_theta` matrix of conditional means.
