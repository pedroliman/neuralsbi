# Sequential NPE with truncated-prior proposals (TSNPE)

Multi-round NPE targeting a single observation `x_obs`. Single-round
[`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) spends its
simulation budget across the whole prior; when only one observation
matters, most of those simulations land in regions the posterior never
visits. `npe_sequential()` implements truncated sequential NPE (TSNPE,
Deistler et al. 2022): after each round the prior is truncated to the
highest-probability region of the current posterior estimate, and the
next round's parameters are drawn from that truncated prior. Because
every proposal is proportional to the prior on its support, the standard
NPE loss stays valid – no importance-weight or atomic correction is
needed, which is what makes TSNPE the simplest correct sequential
scheme.

## Usage

``` r
npe_sequential(
  prior,
  simulator,
  x_obs,
  n_rounds = 2L,
  n_simulations = 1000L,
  sim_args = list(),
  density_estimator = c("maf", "mdn", "nsf", "linear_gaussian"),
  epsilon = 1e-04,
  n_truncation_samples = 5000L,
  max_proposal_batches = 200L,
  seed = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- prior:

  An `nsbi_prior` (see
  [`prior_uniform()`](https://neuralsbi.pedrodelima.com/reference/prior_uniform.md),
  [`prior_normal()`](https://neuralsbi.pedrodelima.com/reference/prior_normal.md)).

- simulator:

  A function called once per parameter set, returning one simulated
  observation; see
  [nsbi_simulator](https://neuralsbi.pedrodelima.com/reference/nsbi_simulator.md).

- x_obs:

  The observation to target, as one numeric vector, one-row matrix or
  one-row data frame whose width matches the simulator's output.
  Sequential inference concentrates simulations around the posterior for
  this observation.

- n_rounds:

  Number of rounds, at least 1. Round 1 is ordinary single-round NPE.

- n_simulations:

  Simulation budget per round; either a scalar or a vector of length
  `n_rounds`.

- sim_args:

  Named list of extra arguments passed to every simulator call; see
  [nsbi_simulator](https://neuralsbi.pedrodelima.com/reference/nsbi_simulator.md).

- density_estimator:

  Passed to
  [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) each
  round.

- epsilon:

  Mass cut for the truncation: the proposal region is the `1 - epsilon`
  highest-probability region of the current posterior. Must be a single
  number strictly between 0 and 1.

- n_truncation_samples:

  Posterior draws used to locate the truncation threshold each round.

- max_proposal_batches:

  Cap on rejection-sampling batches per round.

- seed:

  Optional integer seed for reproducibility.

- verbose:

  Print per-round progress.

- ...:

  Passed to
  [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md)
  (estimator and training settings). Checked before round 1 simulates,
  the same as the arguments above.

## Value

An object of class `c("nsbi_snpe", "nsbi_npe")` with a `rounds` field
recording per-round budgets, acceptance rates, and thresholds.

## Details

The rounds accumulate: each round's estimator is trained on all
simulations so far. The final fit is returned as an `nsbi_npe` (subclass
`nsbi_snpe`) and works with
[`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md),
[`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md) and
the diagnostics, but unlike single-round NPE it is *not* amortized: it
is only trustworthy at (or very near) `x_obs`.

Proposal draws are obtained by rejection: prior candidates are kept when
their posterior log-density clears the `epsilon`-quantile threshold of
the current posterior's own draws. If the posterior is much narrower
than the prior the acceptance rate falls; the round then stops after
`max_proposal_batches` batches and continues with the draws it has, with
a warning.

## References

Deistler, Goncalves & Macke (2022), "Truncated proposals for scalable
and hassle-free simulation-based inference", NeurIPS.
[doi:10.48550/arXiv.2210.04815](https://doi.org/10.48550/arXiv.2210.04815)

## Examples

``` r
prior <- prior_normal(mean = c(mu = 0, nu = 0), sd = 1)
simulator <- function(mu, nu) c(mu, nu) + rnorm(2, sd = 0.3)
fit <- npe_sequential(prior, simulator, x_obs = c(0.5, -0.5),
                      n_rounds = 2, n_simulations = 1000,
                      density_estimator = "linear_gaussian")
post <- posterior(fit, x_obs = c(0.5, -0.5))
draws <- sample(post, 1000)
```
