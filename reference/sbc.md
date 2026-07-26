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
  sim_args = list(),
  seed = NULL
)
```

## Arguments

- fit:

  An `nsbi_npe` fit (amortized posterior).

- simulator:

  The simulator used for inference; called once per trial (see
  [nsbi_simulator](https://neuralsbi.pedrodelima.com/reference/nsbi_simulator.md)).

- prior:

  The prior used for inference (defaults to `fit$prior`).

- n_sbc:

  Number of SBC trials (fresh (theta, x) pairs).

- n_posterior_samples:

  Posterior draws per trial (rank resolution).

- sim_args:

  Named list of extra arguments passed to every simulator call; see
  [nsbi_simulator](https://neuralsbi.pedrodelima.com/reference/nsbi_simulator.md).

- seed:

  Optional seed.

## Value

An object of class `nsbi_sbc` with the rank matrix and a per-parameter
uniformity test.

## Details

A trial whose simulation returns non-finite output is dropped, which
lowers the effective `n_sbc`.

The `n_sbc` simulations run across `future` workers when a plan is set
(see
[nsbi_parallel](https://neuralsbi.pedrodelima.com/reference/nsbi_parallel.md));
the ranking loop that follows calls the trained network and always runs
locally.
