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

## Details

When `d` carries an `n_evals` attribute – `slice_sample_nle()` sets one,
since it is the cost of the run the adapted slice width is trying to
keep down (see
[nsbi_mcmc](https://neuralsbi.pedrodelima.com/reference/nsbi_mcmc.md)) –
it is appended to the summary. A `"stan"` sampler's diagnostics have no
such concept and print without it.
