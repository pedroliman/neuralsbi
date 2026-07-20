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
  theta = NULL,
  x = NULL,
  density_estimator = c("mdn", "maf", "nsf", "linear_gaussian"),
  n_components = 5L,
  n_transforms = 5L,
  hidden = c(50L, 50L),
  max_epochs = 500L,
  batch_size = 100L,
  lr = 5e-04,
  validation_fraction = 0.1,
  patience = 20L,
  n_restarts = 1L,
  clip_grad_norm = 5,
  standardize = TRUE,
  seed = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- prior:

  An `nsbi_prior` (see
  [`prior_uniform()`](https://pedroliman.github.io/neural.sbi/reference/prior_uniform.md),
  [`prior_normal()`](https://pedroliman.github.io/neural.sbi/reference/prior_normal.md)).

- simulator:

  A function mapping an `n x dim` matrix of parameters to an `n x d`
  matrix of simulated data. Ignored if `theta` and `x` are given.

- n_simulations:

  Number of prior draws to simulate when `simulator` is used and
  `theta`/`x` are not supplied.

- theta, x:

  Optional pre-computed simulations. If supplied, `simulator` and
  `n_simulations` are ignored.

- density_estimator:

  One of `"mdn"` (neural Mixture Density Network, needs `torch`),
  `"maf"` (Masked Autoregressive Flow, needs `torch`), `"nsf"` (Neural
  Spline Flow, needs `torch`), or `"linear_gaussian"` (closed-form
  baseline, no `torch`), or a function `function(theta, x)` returning a
  fitted estimator.

- n_components, hidden:

  MDN settings: number of mixture components and a vector of
  hidden-layer widths.

- n_transforms:

  MAF/NSF setting: number of stacked autoregressive transforms.

- max_epochs, batch_size, lr, validation_fraction, patience:

  Neural training controls (Adam optimizer, early stopping on validation
  loss).

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

- seed:

  Optional integer seed for reproducibility.

- verbose:

  Print training progress.

- ...:

  Passed to the density estimator.

## Value

An object of class `nsbi_npe`. Turn it into a usable posterior with
[`posterior()`](https://pedroliman.github.io/neural.sbi/reference/posterior.md),
or sample directly with
[`sample()`](https://pedroliman.github.io/neural.sbi/reference/sample.md).

## Examples

``` r
prior <- prior_uniform(c(-2, -2, -2), c(2, 2, 2))
simulator <- function(theta) theta + 1 + matrix(rnorm(length(theta), sd = 0.1),
                                                 nrow = nrow(theta))
fit <- npe(prior, simulator, n_simulations = 2000,
           density_estimator = "linear_gaussian")
post <- posterior(fit, x_obs = c(0.8, 0.6, 0.4))
draws <- sample(post, 1000)
```
