# Resolve the observation a posterior conditions on

[`resolve_x()`](https://neuralsbi.pedrodelima.com/reference/resolve_x.md)
and
[`resolve_x_iid()`](https://neuralsbi.pedrodelima.com/reference/resolve_x_iid.md)
are both thin wrappers around this: an NPE fit maps one observation to
one posterior, so
[`resolve_x()`](https://neuralsbi.pedrodelima.com/reference/resolve_x.md)
calls this with `first_row = TRUE` and keeps only row 1; an NLE fit's
log-likelihood sums over rows as independent observations of the same
parameter, so
[`resolve_x_iid()`](https://neuralsbi.pedrodelima.com/reference/resolve_x_iid.md)
calls this with `first_row = FALSE` and keeps every row. The same
`x_obs` therefore means "the first observation" to
[`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) and "200
observations" to
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md), and a
user moving a working call from one to the other would otherwise get a
posterior conditioned on a single data point with nothing said about it.
Truncating to row 1 warns rather than fails: taking row 1 of a
simulation matrix is a reasonable thing to ask for.

## Usage

``` r
resolve_obs(post, x, first_row, arg = if (first_row) "x" else "obs")
```

## Arguments

- post:

  An `nsbi_posterior` object.

- x:

  Observation to condition on, or `NULL` to use the posterior's `x_obs`.

- first_row:

  `TRUE` to keep row 1 only, with a warning when there was more than one
  row (the
  [`resolve_x()`](https://neuralsbi.pedrodelima.com/reference/resolve_x.md)
  behavior); `FALSE` to keep every row (the
  [`resolve_x_iid()`](https://neuralsbi.pedrodelima.com/reference/resolve_x_iid.md)
  behavior).

- arg:

  Name the caller's argument goes by, for
  [`check_numeric()`](https://neuralsbi.pedrodelima.com/reference/check_numeric.md)'s
  and
  [`check_finite()`](https://neuralsbi.pedrodelima.com/reference/check_finite.md)'s
  error messages. Defaults to `"x"` when `first_row` and `"obs"`
  otherwise, matching
  [`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md)'s
  and
  [`log_prob()`](https://neuralsbi.pedrodelima.com/reference/log_prob.md)'s
  parameter names; the "no observation supplied" message is fixed to the
  same names and does not vary with `arg`.

## Value

A matrix with `dim_x` columns: one row if `first_row`, every row
otherwise.

## Details

A non-finite entry is a different matter and stops here regardless of
`first_row`. An `NA` passes through `apply_standardizer()` into
`de_sample()` (NPE) or into
[`mcmc_init()`](https://neuralsbi.pedrodelima.com/reference/mcmc_init.md)
(NLE) and comes back as unusable output with nothing to say about the
observation – for NPE the complaint lands in
[`stats::quantile()`](https://rdrr.io/r/stats/quantile.html) inside
[`summary()`](https://rdrr.io/r/base/summary.html), for NLE in
[`mcmc_init()`](https://neuralsbi.pedrodelima.com/reference/mcmc_init.md)
complaining about initialization. There is nothing sensible to condition
on either way, so this errors rather than warns.
