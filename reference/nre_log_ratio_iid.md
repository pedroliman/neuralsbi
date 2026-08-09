# The `n_theta x n_obs` matrix of log ratios

[`de_log_lik_iid()`](https://neuralsbi.pedrodelima.com/reference/de_log_lik_iid.md)'s
counterpart for a ratio estimator, and like
[`nre_iid_evaluator()`](https://neuralsbi.pedrodelima.com/reference/nre_iid_evaluator.md)
a plain function rather than a generic: there is no fast path for any
classifier to specialize.

## Usage

``` r
nre_log_ratio_iid(de, x, theta, max_batch = 1e+05)
```

## Arguments

- de:

  A fitted density estimator.

- x:

  Standardized observations, `n_obs x dim_x`.

- theta:

  Standardized parameters, `n_theta x dim_theta`.

- max_batch:

  Largest number of pairs evaluated at once.

## Value

An `n_theta x n_obs` matrix of log ratios.
