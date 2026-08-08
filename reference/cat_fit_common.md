# Print the fit-summary block shared by `print.nsbi_npe()`, `print.nsbi_nle()` and `print.nsbi_snpe()`

Parameter names, outcome names, the embedding line (when the fit used
one), the simulation count and any drops, the best validation loss, and
the dead-network warning read the same fields regardless of which
factorization was learned. `save_fn_name` is the one line that
legitimately differs between an NPE and an NLE fit, so it is the one
argument callers must supply.

## Usage

``` r
cat_fit_common(x, save_fn_name, data_suffix = "")
```

## Arguments

- x:

  An `nsbi_npe`- or `nsbi_nle`-family fit.

- save_fn_name:

  Name of the save function to point to in the dead-network warning,
  e.g. `"save_npe"`; the matching `load_*()` name is derived from it.

- data_suffix:

  Text appended to the "data (dim)" line before its newline, e.g.
  `" per observation"` for
  [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md).
