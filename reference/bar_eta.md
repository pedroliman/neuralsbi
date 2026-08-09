# ETA in seconds from a steps/sec rate

ETA in seconds from a steps/sec rate

## Usage

``` r
bar_eta(rate, pos, total)
```

## Arguments

- rate:

  Steps per second, as returned by
  [`bar_rate()`](https://neuralsbi.pedrodelima.com/reference/bar_rate.md).

- pos, total:

  Current position and (possibly projected) total.

## Value

Seconds remaining, or `NA_real_` if `rate` is not usable.
