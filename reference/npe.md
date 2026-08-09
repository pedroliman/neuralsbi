# Neural Posterior Estimation (NPE)

`npe()` is the main entry point. Given a prior and either a simulator
(which it will call) or a set of pre-computed simulations `(theta, x)`,
it trains a conditional density estimator whose output directly
approximates the posterior \\p(\theta \mid x)\\. This is single-round,
*amortized* NPE: after training once, you can condition on any
observation without re-simulating.

## Usage

``` r
npe(
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
  embedding_net = NULL,
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
  default, matching Python `sbi`), `"mdn"` (neural Mixture Density
  Network, needs `torch`), `"nsf"` (Neural Spline Flow, needs `torch`),
  or `"linear_gaussian"` (closed-form baseline, no `torch`), or a
  function `function(theta, x)` returning a fitted estimator.

- n_components, hidden:

  MDN settings: number of mixture components (default 10, as in `sbi`)
  and a vector of hidden-layer widths.

- n_transforms:

  MAF/NSF setting: number of stacked autoregressive transforms (default
  5, as in `sbi`).

- n_bins, tail_bound:

  NSF settings: number of spline bins per transform (at least 2, since
  the spline needs an interior derivative to fit) and the half-width of
  the interval the spline acts on (outside it the transform is the
  identity).

- embedding_net:

  Optional summary network built with
  [`embedding_mlp()`](https://neuralsbi.pedrodelima.com/reference/embedding_mlp.md).
  When supplied, the neural estimators condition on the learned features
  \\f\_\psi(x)\\ instead of the raw data, training the embedding
  jointly. Ignored (with a warning) by `"linear_gaussian"`.

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

An object of class `nsbi_npe`. Turn it into a usable posterior with
[`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md),
or sample directly with
[`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md).
Save it to disk with
[`save_npe()`](https://neuralsbi.pedrodelima.com/reference/save_npe.md):
a torch-backed fit does not survive
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html).

## Parallel simulation and progress

The simulator runs sequentially unless you declare a future plan –
`library(future); plan(multisession)` – in which case the simulations
are spread across workers. Simulation and training both report progress
with an ETA. See
[nsbi_parallel](https://neuralsbi.pedrodelima.com/reference/nsbi_parallel.md)
and
[nsbi_progress](https://neuralsbi.pedrodelima.com/reference/nsbi_progress.md).

## Examples

``` r
prior <- prior_uniform(c(mu = -2, nu = -2), c(mu = 2, nu = 2))
simulator <- function(mu, nu) c(a = mu + rnorm(1, sd = 0.1),
                                b = nu + rnorm(1, sd = 0.1))
fit <- npe(prior, simulator, n_simulations = 2000,
           density_estimator = "linear_gaussian")
post <- posterior(fit, x_obs = c(0.8, 0.6))
draws <- sample(post, 1000)
```
