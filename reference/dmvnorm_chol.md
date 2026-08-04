# Multivariate normal log density using a precomputed upper-Cholesky factor (`R` such that `Sigma = t(R) %*% R`, i.e. `chol(Sigma)`).

Multivariate normal log density using a precomputed upper-Cholesky
factor (`R` such that `Sigma = t(R) %*% R`, i.e. `chol(Sigma)`).

Always returns the log density: every call site wants
[`log_prob()`](https://neuralsbi.pedrodelima.com/reference/log_prob.md)'s
contract, none of [`dnorm()`](https://rdrr.io/r/stats/Normal.html)'s
`log = FALSE`, so there is no `log` argument to forget to set.

## Usage

``` r
dmvnorm_chol(x, mean, R)
```
