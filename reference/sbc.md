# Simulation-Based Calibration (SBC)

Repeatedly draws a "true" parameter from the prior, simulates data, and
ranks the true parameter within posterior samples conditioned on that
data. If the posterior is well calibrated, the ranks are uniformly
distributed.

## Usage

``` r
sbc(
  fit,
  simulator,
  prior = fit$prior,
  n_sbc = 200L,
  n_posterior_samples = 1000L,
  seed = NULL
)
```

## Arguments

- fit:

  An `nsbi_npe` fit (amortized posterior).

- simulator:

  The simulator used for inference.

- prior:

  The prior used for inference (defaults to `fit$prior`).

- n_sbc:

  Number of SBC trials (fresh (theta, x) pairs).

- n_posterior_samples:

  Posterior draws per trial (rank resolution).

- seed:

  Optional seed.

## Value

An object of class `nsbi_sbc` with the rank matrix and a per-parameter
uniformity test.
