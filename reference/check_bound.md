# Validate a support bound

[`within_support()`](https://neuralsbi.pedrodelima.com/reference/within_support.md)
compares `theta` against `lower`/`upper` with
[`sweep()`](https://rdrr.io/r/base/sweep.html), which recycles a bound
of the wrong length and warns instead of stopping. The support test that
comes back is then wrong, and it decides which posterior draws are
rejected as leakage and what
[`log_prob()`](https://neuralsbi.pedrodelima.com/reference/log_prob.md)
renormalizes by. A length-1 bound is recycled here, once, so everything
downstream sees one bound per parameter.

## Usage

``` r
check_bound(value, arg, d)
```

## Arguments

- value:

  The user's value, or `NULL` for an unbounded side.

- arg:

  Name of the argument.

- d:

  Number of parameters.

## Value

`NULL`, or a double vector of length `d`.
