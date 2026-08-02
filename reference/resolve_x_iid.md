# All rows of the observation, unlike `resolve_x()` which keeps only the first

A non-finite entry is rejected here for the same reason it is in
[`resolve_x()`](https://neuralsbi.pedrodelima.com/reference/resolve_x.md):
the log-likelihood sums over rows, so one `NA` makes every starting
point non-finite and the run fails in
[`mcmc_init()`](https://neuralsbi.pedrodelima.com/reference/mcmc_init.md)
complaining about initialization instead of about the observation.

## Usage

``` r
resolve_x_iid(post, x, arg = "obs")
```

## Arguments

- post:

  An `nsbi_nle_posterior` object.

- x:

  Observation to condition on, or `NULL` to use the posterior's `x_obs`.

- arg:

  Name the caller's argument goes by, for the error message.
  [`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md)
  calls it `obs` and
  [`log_prob()`](https://neuralsbi.pedrodelima.com/reference/log_prob.md)
  calls it `x`.
