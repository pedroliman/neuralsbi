# Neural Likelihood Estimation (NLE)

`nle()` trains a conditional density estimator on the *other*
factorization of the joint: instead of the posterior \\p(\theta \mid
x)\\ that [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md)
learns, it learns a surrogate likelihood \\q\_\phi(x \mid \theta)\\. The
posterior then follows from Bayes' rule, \\p(\theta \mid x) \propto
q\_\phi(x \mid \theta)\\p(\theta)\\, and is sampled with MCMC by
[`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md).

## Usage

``` r
nle(
  prior,
  simulator = NULL,
  n_simulations = 1000,
  sim_args = list(),
  theta = NULL,
  x = NULL,
  density_estimator = c("maf", "mdn", "nsf", "linear_gaussian"),
  n_components = 10L,
  n_transforms = 5L,
  hidden = c(50L, 50L),
  n_bins = 10L,
  tail_bound = 3,
  max_epochs = 2000L,
  batch_size = 200L,
  lr = 5e-04,
  validation_fraction = 0.1,
  patience = 20L,
  n_restarts = 1L,
  clip_grad_norm = 5,
  standardize = TRUE,
  device = "cpu",
  seed = NULL,
  verbose = FALSE
)
```

## Arguments

- prior:

  An `nsbi_prior` (see
  [`prior_uniform()`](https://neuralsbi.pedrodelima.com/reference/prior_uniform.md),
  [`prior_normal()`](https://neuralsbi.pedrodelima.com/reference/prior_normal.md)).

- simulator:

  A function called once per parameter set, returning one simulated
  observation: a numeric vector, a scalar, or a one-row matrix or data
  frame. Parameters arrive either as named arguments (when the prior's
  names match the simulator's formals) or as one named vector. Names on
  the output become the outcome names used in plots. See
  [nsbi_simulator](https://neuralsbi.pedrodelima.com/reference/nsbi_simulator.md).
  Ignored if `theta` and `x` are given.

- n_simulations:

  Number of prior draws to simulate when `simulator` is used and
  `theta`/`x` are not supplied.

- sim_args:

  Named list of extra arguments passed to every simulator call: observed
  data, a time grid, a fixed population size, solver settings. See
  [nsbi_simulator](https://neuralsbi.pedrodelima.com/reference/nsbi_simulator.md).

- theta, x:

  Optional pre-computed simulations. If supplied, `simulator` and
  `n_simulations` are ignored. Column names on `theta` (or names on
  `prior`'s `mean`/`low`) and on `x` are carried through to posterior
  samples, SBC results, and their plots.

- density_estimator:

  One of `"maf"` (Masked Autoregressive Flow, needs `torch`; the
  default), `"mdn"` (Mixture Density Network, needs `torch`), `"nsf"`
  (Neural Spline Flow, needs `torch`), or `"linear_gaussian"`
  (closed-form baseline, no `torch`), or a function `function(theta, x)`
  returning a fitted estimator. Note the estimator sees the roles
  swapped: its target is `x` and it conditions on `theta`.

- n_components, hidden:

  MDN settings: number of mixture components and a vector of
  hidden-layer widths.

- n_transforms:

  MAF/NSF setting: number of stacked autoregressive transforms.

- n_bins, tail_bound:

  NSF settings.

- max_epochs, batch_size, lr, validation_fraction, patience:

  Neural training controls (Adam optimizer, early stopping on validation
  loss). The defaults (`batch_size = 200`, `lr = 5e-4`,
  `validation_fraction = 0.1`, `patience = 20`) match Python `sbi`;
  `max_epochs` is a high guard cap that early stopping normally reaches
  first.

- n_restarts:

  Train this many independently initialized networks and keep the one
  with the best validation loss (guards against bad initializations and
  MDN mode collapse).

- clip_grad_norm:

  Maximum gradient norm during training (`Inf` disables clipping). The
  learning rate also decays 2x after 10 epochs without validation
  improvement.

- standardize:

  Whether to z-score `theta` and `x` before training (strongly
  recommended; default `TRUE`).

- device:

  Where to train the neural estimator: `"cpu"` (the default), `"cuda"`,
  `"mps"`, or `"gpu"`/`"auto"` to resolve CUDA -\> MPS -\> CPU
  (mirroring Python `sbi`'s `"gpu"`). CPU is the default on purpose –
  matching `sbi`, not auto-selecting a GPU – and `"cuda"`/`"mps"` error
  if the requested device is not actually available rather than falling
  back silently, so a typo or a missing driver is not mistaken for a
  slow CPU run. Only `"gpu"`/`"auto"` falls back to CPU without
  complaint, since it never named a specific device. Ignored (with no
  error) by `"linear_gaussian"`, which has no GPU concept. For the small
  networks typical of SBI models (the SIR example, say),
  `"mps"`/`"cuda"` can be *slower* than CPU – per-kernel launch overhead
  dominates until the net and batch are large – so try CPU first and
  switch only if profiling shows a win.

- seed:

  Optional integer seed for reproducibility.

- verbose:

  Print training progress.

## Value

An object of class `nsbi_nle`. Evaluate the surrogate likelihood with
[`log_lik()`](https://neuralsbi.pedrodelima.com/reference/log_lik.md),
turn it into a posterior with
[`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md),
or export it to Stan with
[`stan_code()`](https://neuralsbi.pedrodelima.com/reference/stan_export.md).

## When NLE beats NPE

An NPE fit learns the posterior for one fixed data dimension, chosen at
training time. If the observation is \\n\\ independent trials from the
same parameter, NPE must be retrained for every \\n\\, or handed summary
statistics that throw information away. NLE learns the density of a
*single* trial, so the log-likelihood of \\n\\ trials is a sum of \\n\\
evaluations:

\$\$\log p(x_1, \ldots, x_n \mid \theta) = \sum\_{i=1}^{n} \log
q\_\phi(x_i \mid \theta).\$\$

Train once, then condition on 50 trials or 5000 without touching the
network. The learned likelihood is also a plain differentiable function
of \\\theta\\, so it can be embedded in a larger model written by hand –
see
[`stan_code()`](https://neuralsbi.pedrodelima.com/reference/stan_export.md).

The trade-off is real and worth stating: posterior draws now cost an
MCMC run rather than a forward pass, and for a single fixed observation
with high-dimensional data NPE is usually the better choice.

## References

Papamakarios, G., Sterratt, D. and Murray, I. (2019). Sequential Neural
Likelihood: Fast Likelihood-free Inference with Autoregressive Flows.
*AISTATS*.
[doi:10.48550/arXiv.1805.07226](https://doi.org/10.48550/arXiv.1805.07226)

## See also

[`log_lik()`](https://neuralsbi.pedrodelima.com/reference/log_lik.md)
and
[`likelihood_fn()`](https://neuralsbi.pedrodelima.com/reference/likelihood_fn.md)
to evaluate the surrogate,
[`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md)
to sample it,
[`stan_code()`](https://neuralsbi.pedrodelima.com/reference/stan_export.md)
to hand it to Stan.

## Examples

``` r
# One noisy measurement per simulator call; the observation is 200 of them.
prior <- prior_uniform(c(mu = -3, log_sigma = -1),
                       c(mu =  3, log_sigma =  1))
simulator <- function(mu, log_sigma) c(y = rnorm(1, mu, exp(log_sigma)))

fit <- nle(prior, simulator, n_simulations = 2000,
           density_estimator = "linear_gaussian")

x_obs <- matrix(rnorm(200, mean = 1, sd = 0.5), ncol = 1)
log_lik(fit, theta = c(1, log(0.5)), x = x_obs)
#> [1] -255.1691
```
