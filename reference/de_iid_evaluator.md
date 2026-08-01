# A summed i.i.d. log-likelihood with the observation held fixed

[`log_lik()`](https://neuralsbi.pedrodelima.com/reference/log_lik.md)
and every MCMC step ask the same question over and over: the summed
log-density of one fixed set of observations under a `theta` that
changes. `de_iid_evaluator()` returns a closure over the observations,
so whatever an estimator can settle once settles when the closure is
built rather than on every call. For the MDN that is coercing the
observations to a tensor, which at a few thousand rows is not a rounding
error next to the forward pass.

## Usage

``` r
de_iid_evaluator(de, x, max_batch = 1e+05)
```

## Arguments

- de:

  A fitted density estimator.

- x:

  Standardized observations, `n_obs x dim_x`.

- max_batch:

  Largest number of pairs evaluated at once.

## Value

`function(theta)` giving one summed log-density per row of `theta`, in
standardized space.

## Details

Reducing inside the closure matters as much as the hoisting. The
`n_theta x n_obs` matrix is the largest object in the loop and none of
it is wanted, so the sum happens where the log-densities are produced
and only `n_theta` numbers ever cross back into R.

Estimators need a method here only if they can beat the default, which
is
[`de_log_lik_iid()`](https://neuralsbi.pedrodelima.com/reference/de_log_lik_iid.md)
with its row sums taken block by block.
