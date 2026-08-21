# Draw from one marginal by inverting its CDF between the bounds

Inverse-CDF sampling covers the truncated and untruncated cases with the
same line, and it is exact: no rejection loop that stalls when the
bounds cut into a tail.

## Usage

``` r
marginal_sample(m, n)
```

## Arguments

- m:

  A marginal from
  [`new_marginal()`](https://neuralsbi.pedrodelima.com/reference/new_marginal.md).

- n:

  Number of draws.
