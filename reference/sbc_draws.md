# Shared simulate-and-drop preamble for sbc() and tarp()

Both diagnostics draw "true" parameters from `prior`, simulate data for
them, and drop trials whose simulation failed. The prior-width check
lives here rather than in each caller, so
[`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md) and
[`tarp()`](https://neuralsbi.pedrodelima.com/reference/tarp.md) get it
in one place: `prior` exists to be overridden (see their docs), and a
prior of the wrong width would otherwise reach
[`sweep()`](https://rdrr.io/r/base/sweep.html)/z-scoring downstream and
silently recycle against the fit's width.

## Usage

``` r
sbc_draws(fit, simulator, prior, n, sim_args, what)
```

## Arguments

- fit:

  An `nsbi_npe` or `nsbi_nle` fit.

- simulator:

  The simulator used for inference.

- prior:

  The prior to draw true parameters from.

- n:

  Number of trials to draw.

- sim_args:

  Named list of extra simulator arguments.

- what:

  Label for
  [`drop_failed_sims()`](https://neuralsbi.pedrodelima.com/reference/drop_failed_sims.md)'s
  warning (e.g. `"SBC trials"`).

## Value

A list with `theta` and `x` (the surviving trials), `n` (their count,
i.e. `nrow(theta)` after dropping), and `n_dropped`.
