# The single observation an NPE posterior conditions on

Thin wrapper around
[`resolve_obs()`](https://neuralsbi.pedrodelima.com/reference/resolve_obs.md)
with `first_row = TRUE`; see there for the rationale shared with
[`resolve_x_iid()`](https://neuralsbi.pedrodelima.com/reference/resolve_x_iid.md).

## Usage

``` r
resolve_x(post, x)
```

## Arguments

- post:

  An `nsbi_posterior` object.

- x:

  Observation to condition on, or `NULL` to use the posterior's `x_obs`.

## Value

A one-row matrix with `dim_x` columns.
