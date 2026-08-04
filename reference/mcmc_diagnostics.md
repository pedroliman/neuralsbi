# Split-Rhat and bulk effective sample size

The standard rank-free versions from Vehtari et al. (2021), computed on
the split chains. Implemented here rather than taken from posterior to
keep the dependency surface where it is; the test suite cross-checks
against posterior when that package happens to be installed.

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
