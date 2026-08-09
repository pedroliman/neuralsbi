# Recent-window step rate for the built-in bar's ETA

`elapsed / pos` – the lifetime average since the bar started – is what
the old ETA extrapolated from, and it is a poor estimate of the
*current* speed: it is dragged down for the whole run by the first
step's one-time setup cost (device/net initialization, first CUDA/MPS
kernel compilation, R JIT warmup), and across training restarts it
blends restarts that may cost very different amounts per epoch.
`times`/`positions` are a short rolling window of recent calls (oldest
first), so the rate reflects what the last few steps actually cost, not
the whole run's history. Falls back to the lifetime average when the
window is too short to trust (fewer than two samples, or no
time/progress elapsed within it yet).

## Usage

``` r
bar_rate(times, positions, elapsed, pos)
```

## Arguments

- times, positions:

  Parallel vectors of recent
  [`Sys.time()`](https://rdrr.io/r/base/Sys.time.html) timestamps (as
  numeric seconds) and the `pos` reported at each.

- elapsed, pos:

  Lifetime elapsed seconds and current position, used only as the
  fallback.

## Value

Steps per second, or `NA_real_` if no rate can be estimated yet.
