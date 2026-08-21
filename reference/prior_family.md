# The distribution families a prior can be built from

One entry per family: the parameter names in Stan's order, the R
density/CDF/quantile trio behind them, and the natural support.
Everything else in this file is generic over that entry, which is why
registering a family is all it takes for
[`prior_truncated()`](https://neuralsbi.pedrodelima.com/reference/prior_truncated.md)
to renormalize it and for
[`stan_code()`](https://neuralsbi.pedrodelima.com/reference/stan_export.md)
to write it out.

## Usage

``` r
prior_family(family)
```

## Arguments

- family:

  Family name, or anything else to get `NULL` back.

## Value

A list, or `NULL` for an unregistered name.

## Details

The `d`, `p` and `q` functions are always called with the value
positionally and the distribution parameters by name, since R spells the
first argument `x`, `q` and `p` in turn while the parameters keep their
names throughout.
