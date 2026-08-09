# Minimal dependency-free progress bar with an ETA

Draws to [`stderr()`](https://rdrr.io/r/base/showConnections.html),
redrawing in place. Updates are throttled so a fast inner loop does not
spend its time formatting text. The ETA is extrapolated from a short
rolling window of recent calls (see
[`bar_rate()`](https://neuralsbi.pedrodelima.com/reference/bar_rate.md)),
not the lifetime average, so it tracks the current speed rather than a
blend of setup overhead and past restarts.

## Usage

``` r
builtin_bar(label = NULL)
```
