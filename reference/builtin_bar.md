# Minimal dependency-free progress bar with an ETA

Draws to [`stderr()`](https://rdrr.io/r/base/showConnections.html),
redrawing in place. Updates are throttled so a fast inner loop does not
spend its time formatting text.

## Usage

``` r
builtin_bar(label = NULL)
```
