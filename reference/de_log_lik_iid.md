# Log-density of many observations under many parameter values

The cross product `n_theta x n_obs`, in standardized space. Only
`sum_iid = FALSE` wants that matrix; everything else –
[`log_lik()`](https://neuralsbi.pedrodelima.com/reference/log_lik.md)'s
default and every MCMC step – wants its row sums, and gets them from
[`de_iid_evaluator()`](https://neuralsbi.pedrodelima.com/reference/de_iid_evaluator.md)
without ever building it.

## Usage

``` r
de_log_lik_iid(de, x, theta, max_batch = 1e+05)
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

An `n_theta x n_obs` matrix of log-densities.

## Details

The default expands the cross product and makes one batched call, which
is what a flow needs: its transforms depend on the observation as well
as the parameter, so there is nothing to reuse between observations.
Estimators whose conditional distribution depends on the parameter
*alone* – the MDN and the linear-Gaussian baseline – override this and
compute that distribution once per parameter, which turns the i.i.d. sum
from `n_theta * n_obs` network passes into `n_theta` of them. With a few
thousand observations that is the difference between usable and not.
