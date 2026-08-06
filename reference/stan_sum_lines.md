# The i.i.d.-sum entry point shared by every generated `_sum_lpdf`

Standardizes `theta`, declares `x`'s center/scale, and accumulates
`body` (an expression for one observation's log density, in terms of
`x[n]`) over `rows(x)` before applying the jacobian once at the end.
`precompute` is the one place estimators differ: `linear_gaussian` and
the MDN can build their conditional distribution once, outside the loop,
because it depends on `theta` alone; MAF has nothing to hoist, so its
caller leaves this at the default. Mirrors
[`de_log_lik_iid()`](https://neuralsbi.pedrodelima.com/reference/de_log_lik_iid.md)
(R/likelihood.R), which sums the same per-observation log density on the
R side; the two have to agree.

## Usage

``` r
stan_sum_lines(fit, P, body, precompute = "")
```
