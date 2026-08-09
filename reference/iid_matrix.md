# The `n_theta x n_obs` cross product, block by block

The body of
[`de_log_lik_iid()`](https://neuralsbi.pedrodelima.com/reference/de_log_lik_iid.md)'s
default method, with the per-pair scorer left open so
[`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md) can reuse
it. `score(de, x_rows, theta_rows)` is
[de_log_prob()](https://neuralsbi.pedrodelima.com/reference/density_estimator.md)
for a density estimator and
[`nre_score()`](https://neuralsbi.pedrodelima.com/reference/nre_score.md)
for a ratio estimator.

## Usage

``` r
iid_matrix(de, x, theta, max_batch, score)
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

- score:

  The per-pair scorer, `function(de, x, theta)`.
