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
  sim_args = list(),
  seed = NULL,
  ...
)
```

## Arguments

- fit:

  An `nsbi_npe` fit from
  [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md), or an
  `nsbi_nle` fit from
  [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md). With an
  NLE fit every trial is a separate MCMC run, so start with a small
  `n_tarp` and raise it once the cost is known.

- simulator:

  The simulator used for inference; called once per trial (see
  [nsbi_simulator](https://neuralsbi.pedrodelima.com/reference/nsbi_simulator.md)).

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

- sim_args:

  Named list of extra arguments passed to every simulator call; see
  [nsbi_simulator](https://neuralsbi.pedrodelima.com/reference/nsbi_simulator.md).

- seed:

  Optional seed.

- ...:

  Passed to
  [`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md),
  which is how the MCMC controls (`n_chains`, `warmup`, `thin`,
  `sampler`) reach an NLE fit.

## Value

An object of class `nsbi_tarp` with the per-trial coverage values and
the ECP curve. Plot it with
[`plot_tarp()`](https://neuralsbi.pedrodelima.com/reference/plot_tarp.md).

## Details

Unlike [`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md),
which ranks each parameter marginally, TARP is a *joint* test: it can
detect posteriors whose marginals are calibrated but whose correlation
structure is wrong. Distances are computed after z-scoring each
parameter (using the spread of the true draws), so parameters on
different scales contribute comparably.

A trial whose simulation returns non-finite output is dropped, which
lowers the effective `n_tarp`. A trial whose posterior comes back with
fewer than `n_posterior_samples` draws is an error, for the same reason
as in [`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md): a
trial scored on a different number of draws is not comparable to the
rest.

## References

Lemos, Coogan, Hezaveh & Perreault-Levasseur (2023), "Sampling-based
accuracy testing of posterior estimators for general inference", ICML.
[doi:10.48550/arXiv.2302.03026](https://doi.org/10.48550/arXiv.2302.03026)
