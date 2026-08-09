# The body [`log_lik()`](https://neuralsbi.pedrodelima.com/reference/log_lik.md) and [`log_ratio()`](https://neuralsbi.pedrodelima.com/reference/log_ratio.md) share

Both take a `theta` and an `x` whose rows are independent observations,
check them, standardize them, and either return the `n_theta x n_obs`
matrix or its row sums. Which functions do the scoring, and whether the
Jacobian comes with them, is
[`surrogate_ops()`](https://neuralsbi.pedrodelima.com/reference/surrogate_ops.md)'s
business.

## Usage

``` r
surrogate_score(fit, theta, x, sum_iid, max_batch)
```

## Arguments

- fit:

  An `nsbi_nle` or `nsbi_nre` fit.

- theta, x, sum_iid, max_batch:

  As in
  [`log_lik()`](https://neuralsbi.pedrodelima.com/reference/log_lik.md).
