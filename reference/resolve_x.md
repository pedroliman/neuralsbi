# The single observation an NPE posterior conditions on

An NPE fit maps one observation to one posterior, so only the first row
of `x` can be used.
[`resolve_x_iid()`](https://neuralsbi.pedrodelima.com/reference/resolve_x_iid.md),
the NLE counterpart, keeps every row, because there the rows are
independent observations and the log-likelihood sums over them. The same
`x_obs` therefore means "200 observations" to
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) and "the
first observation" to
[`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md), and a
user moving a working call from one to the other would otherwise get a
posterior conditioned on a single data point with nothing said about it.
Warn rather than fail: taking row 1 of a simulation matrix is a
reasonable thing to ask for.

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
