# Shared per-trial posterior-draw loop for sbc() and tarp()

Draws a posterior for each row of `x_all`, insists on the full
`n_posterior_samples` count via
[`diagnostic_draws()`](https://neuralsbi.pedrodelima.com/reference/diagnostic_draws.md),
and hands the result to `f(draws, i)` for whichever per-trial metric the
caller wants:
[`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md) returns a
rank row,
[`tarp()`](https://neuralsbi.pedrodelima.com/reference/tarp.md) a scalar
coverage value.

## Usage

``` r
for_each_trial(fit, x_all, n_posterior_samples, label, f, ...)
```

## Arguments

- fit:

  An `nsbi_npe`, `nsbi_nle` or `nsbi_nre` fit.

- x_all:

  Matrix of simulated data, one row (one trial) per observation.

- n_posterior_samples:

  Posterior draws requested per trial.

- label:

  Progress-bar label.

- f:

  `function(draws, i)`, the per-trial metric.

- ...:

  Passed to
  [`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md).

## Value

A list with `results` (a list of length `nrow(x_all)`, one `f()` return
value per trial) and `n_posterior_samples` (the draw count trials were
actually scored against).

## Details

The denominator a caller bins or normalizes against
([`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md)'s rank
scale, via the p-value bins and
[`expected_coverage()`](https://neuralsbi.pedrodelima.com/reference/expected_coverage.md))
is read back from how many draws a trial actually returned, not trusted
from the `n_posterior_samples` argument – so if a future change ever
lets a short draw through instead of erroring, the diagnostic is still
scored on the scale it was actually drawn on rather than the one it was
asked for.

`finally`, not [`on.exit()`](https://rdrr.io/r/base/on.exit.html): this
block runs inside
[`progressr::with_progress()`](https://progressr.futureverse.org/reference/with_progress.html),
so [`on.exit()`](https://rdrr.io/r/base/on.exit.html) would attach to
the promise's forcing frame and close the bar before the loop starts
(see
[`run_simulator()`](https://neuralsbi.pedrodelima.com/reference/run_simulator.md)
in R/parallel.R).
