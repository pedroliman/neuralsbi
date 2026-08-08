# Print an "N dropped" line shared by fit summaries and calibration diagnostics

A fit knows its simulation budget, so its line reports a drop rate;
[`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md) and
[`tarp()`](https://neuralsbi.pedrodelima.com/reference/tarp.md) only
know how many further trials were lost mid-run, so leaving `total` at
`NULL` switches to that shorter wording.

## Usage

``` r
cat_dropped(n_dropped, total = NULL, what)
```

## Arguments

- n_dropped:

  Number dropped. Nothing is printed when this is `NULL` or `0`.

- total:

  Simulations attempted (dropped plus kept), or `NULL` for the
  diagnostics wording, which has no rate to report.

- what:

  What was dropped, e.g. `"non-finite output"`.
