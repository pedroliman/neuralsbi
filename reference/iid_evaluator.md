# Summed log densities over a fixed observation set, block by block

The body of
[`de_iid_evaluator()`](https://neuralsbi.pedrodelima.com/reference/de_iid_evaluator.md)'s
default method, with the per-pair scorer left open so
[`nre_iid_evaluator()`](https://neuralsbi.pedrodelima.com/reference/nre_iid_evaluator.md)
can reuse it.

## Usage

``` r
iid_evaluator(de, x, max_batch, score)
```

## Arguments

- de:

  A fitted density estimator.

- x:

  Standardized observations, `n_obs x dim_x`.

- max_batch:

  Largest number of pairs evaluated at once.

- score:

  The per-pair scorer, `function(de, x, theta)`.
