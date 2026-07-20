# TARP expected coverage

Tests of Accuracy with Random Points (Lemos et al. 2023). For each trial
a true parameter is drawn from the prior, data are simulated, and
posterior samples are drawn conditioned on those data. Given a random
reference point, the fraction of posterior samples closer to the
reference than the truth is the credibility level of the smallest
distance-based credible region that contains the truth. For a calibrated
posterior these fractions are uniform, so the expected coverage
probability (ECP) at credibility level alpha equals alpha.

## Usage

``` r
tarp(
  fit,
  simulator,
  prior = fit$prior,
  n_tarp = 200L,
  n_posterior_samples = 1000L,
  references = c("uniform", "prior"),
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

- n_tarp:

  Number of TARP trials (fresh (theta, x) pairs).

- n_posterior_samples:

  Posterior draws per trial.

- references:

  How to draw reference points: `"uniform"` (default, uniform over the
  hyper-rectangle spanned by the true parameter draws, as in the paper)
  or `"prior"` (draws from the prior).

- seed:

  Optional seed.

## Value

An object of class `nsbi_tarp` with the per-trial coverage values and
the ECP curve. Plot it with
[`plot_tarp()`](https://pedroliman.github.io/neural.sbi/reference/plot_tarp.md).

## Details

Unlike
[`sbc()`](https://pedroliman.github.io/neural.sbi/reference/sbc.md),
which ranks each parameter marginally, TARP is a *joint* test: it can
detect posteriors whose marginals are calibrated but whose correlation
structure is wrong. Distances are computed after z-scoring each
parameter (using the spread of the true draws), so parameters on
different scales contribute comparably.

## References

Lemos, Coogan, Hezaveh & Perreault-Levasseur (2023), "Sampling-based
accuracy testing of posterior estimators for general inference", ICML.
[doi:10.48550/arXiv.2302.03026](https://doi.org/10.48550/arXiv.2302.03026)
