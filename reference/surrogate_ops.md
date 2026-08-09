# How a surrogate fit is scored

[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) and
[`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md) fits are
interchangeable everywhere downstream of the estimator, and they differ
in exactly three things: which function produces the `n_theta x n_obs`
matrix of scores, which one produces its row sums without building it,
and whether standardizing `x` needs a change-of-variables term. A
density reported in the simulator's units needs one; a ratio does not,
since the Jacobian cancels between its numerator and denominator (see
[`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md)).

## Usage

``` r
surrogate_ops(fit)
```

## Arguments

- fit:

  An `nsbi_nle` or `nsbi_nre` fit.

## Value

A list with `matrix_fn`, `evaluator` and `log_jac`.

## Details

Keeping the three together in one table is what lets
[`surrogate_score()`](https://neuralsbi.pedrodelima.com/reference/surrogate_score.md)
and
[`surrogate_potential()`](https://neuralsbi.pedrodelima.com/reference/surrogate_potential.md)
be plain functions rather than a generic each. It also means the fact
that a ratio has no Jacobian is written down once.
