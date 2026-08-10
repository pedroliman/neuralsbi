# All rows of the observation, unlike `resolve_x()` which keeps only the first

Thin wrapper around
[`resolve_obs()`](https://neuralsbi.pedrodelima.com/reference/resolve_obs.md)
with `first_row = FALSE`; see there for the rationale shared with
[`resolve_x()`](https://neuralsbi.pedrodelima.com/reference/resolve_x.md),
including why a non-finite entry stops here rather than warns.

## Usage

``` r
resolve_x_iid(post, x, arg = "obs")
```

## Arguments

- post:

  An `nsbi_mcmc_posterior` object.

- x:

  Observation to condition on, or `NULL` to use the posterior's `x_obs`.

- arg:

  Name the caller's argument goes by, for the error message.
  [`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md)
  calls it `obs` and
  [`log_prob()`](https://neuralsbi.pedrodelima.com/reference/log_prob.md)
  calls it `x`.
