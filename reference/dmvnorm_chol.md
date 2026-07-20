# Multivariate normal log density using a precomputed upper-Cholesky factor (`R` such that `Sigma = t(R) %*% R`, i.e. `chol(Sigma)`).

Multivariate normal log density using a precomputed upper-Cholesky
factor (`R` such that `Sigma = t(R) %*% R`, i.e. `chol(Sigma)`).

## Usage

``` r
dmvnorm_chol(x, mean, R, log = TRUE)
```
