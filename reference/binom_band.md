# 99% Monte-Carlo binomial band for a calibration curve

At each nominal level, how far empirical coverage over `n` trials can
wander from the diagonal by chance alone, under a `Binomial(n, nominal)`
model of the count that lands inside its interval.
[`plot_coverage()`](https://neuralsbi.pedrodelima.com/reference/plot_coverage.md)
and
[`plot_tarp()`](https://neuralsbi.pedrodelima.com/reference/plot_tarp.md)
shade this as the ribbon behind their curve;
[`plot_sbc()`](https://neuralsbi.pedrodelima.com/reference/plot_sbc.md)
uses the same interval, scaled back up to a count, for its histogram's
dashed reference lines. One formula in one place means all three
figures' bands agree on what "real" departure from calibration looks
like.

## Usage

``` r
binom_band(nominal, n, level = 0.99)
```

## Arguments

- nominal:

  Nominal level(s) the band is evaluated at.

- n:

  Number of trials the empirical coverage was computed over.

- level:

  Two-sided coverage of the band, e.g. `0.99` for the `c(0.005, 0.995)`
  binomial quantiles.

## Value

A data frame with columns `nominal`, `lo` and `hi`, all as fractions of
`n`.
