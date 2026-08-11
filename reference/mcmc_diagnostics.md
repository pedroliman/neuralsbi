# Split-Rhat and bulk effective sample size

The standard rank-normalized versions from Vehtari et al. (2021),
computed on the split chains. `split_rhat()` is
`max(bulk-Rhat, tail-Rhat)`: the classical Gelman-Rubin statistic run
once on the rank-normalized draws (bulk) and once on the rank-normalized
draws folded around the median (tail, which catches chains that agree in
location but disagree in spread – something bulk-Rhat alone can miss).
`bulk_ess()` rank-normalizes the same way. Implemented here rather than
taken from posterior to keep the dependency surface where it is; the
test suite cross-checks against posterior when that package happens to
be installed.

## Usage

``` r
mcmc_diagnostics(chains)
```

## Arguments

- chains:

  A `n_iter x n_chains x dim` array.

## Value

A data frame with one row per parameter: `rhat` and `ess_bulk`.

## References

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Burkner, P.-C.
(2021). Rank-normalization, folding, and localization. *Bayesian
Analysis* 16(2), 667-718.
[doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)
