# Summed log ratio over independent observations, with the observation fixed

[`de_iid_evaluator()`](https://neuralsbi.pedrodelima.com/reference/de_iid_evaluator.md)'s
counterpart for a ratio estimator. There is no fast path to specialize:
the classifier sees `(theta, x_i)` jointly, so every pair costs a
forward pass however the loop is arranged.

## Usage

``` r
nre_iid_evaluator(de, x, max_batch = 1e+05)
```

## Arguments

- de:

  A fitted density estimator.

- x:

  Standardized observations, `n_obs x dim_x`.

- max_batch:

  Largest number of pairs evaluated at once.

## Value

`function(theta)` giving one summed log ratio per row of `theta`.
