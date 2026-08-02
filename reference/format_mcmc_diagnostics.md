# Summarize MCMC diagnostics for printing

`split_rhat()` and `bulk_ess()` return `NA` for a run they cannot score:
too few iterations, one chain, or a coordinate that never moved. Every
parameter can come back `NA` at once, and `max(na.rm = TRUE)` over
all-`NA` is `-Inf` with a warning, which the print method used to show
as if it were a value. A partially scored run is reported from the
parameters that did get a number, with a count of the ones that did not.

## Usage

``` r
format_mcmc_diagnostics(d)
```

## Arguments

- d:

  The data frame from
  [`mcmc_diagnostics()`](https://neuralsbi.pedrodelima.com/reference/mcmc_diagnostics.md).

## Value

One line of text, without a trailing newline.
