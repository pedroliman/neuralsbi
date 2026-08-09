# Neural Ratio Estimation (NRE)

`nre()` learns the third factorization of the joint.
[`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) learns the
posterior and
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) learns the
likelihood; `nre()` learns neither density but their ratio, \\r(\theta,
x) = p(x \mid \theta) / p(x)\\, by training a binary classifier to tell
\\(\theta, x)\\ pairs drawn from the joint apart from pairs whose
parameter came from a different simulation. The posterior follows from
Bayes' rule, \\p(\theta \mid x) \propto r(\theta, x)\\p(\theta)\\, and
is sampled with MCMC by
[`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md).

## Usage

``` r
nre(
  prior,
  simulator = NULL,
  n_simulations = 1000,
  sim_args = list(),
  theta = NULL,
  x = NULL,
  classifier = c("resnet", "mlp", "linear", "logistic"),
  num_atoms = 10L,
  hidden = 50L,
  n_blocks = 2L,
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

- classifier:

  One of `"resnet"` (a residual MLP, needs `torch`; the default,
  matching Python `sbi`), `"mlp"` (a plain MLP, needs `torch`),
  `"linear"` (a single linear layer on the raw inputs, needs `torch`),
  or `"logistic"` (a closed-form logistic regression on quadratic
  features, no `torch`), or a function `function(theta, x)` returning a
  fitted ratio estimator – one whose class has a
  [`de_log_ratio()`](https://neuralsbi.pedrodelima.com/reference/de_log_ratio.md)
  method, which is the only thing the rest of the pipeline asks of it.

- num_atoms:

  Number of parameter values the classifier compares per simulation: one
  true and `num_atoms - 1` contrasts. Clamped to the minibatch size, as
  in `sbi`.

- hidden:

  Width of the classifier's hidden layers. One number, not a per-layer
  vector as in
  [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md): `sbi`'s
  classifiers use a single width throughout, and `n_blocks` is what sets
  the depth.

- n_blocks:

  Depth of the classifier: residual blocks for `"resnet"`, hidden layers
  for `"mlp"`. Ignored by `"linear"` and `"logistic"`. `sbi` fixes its
  MLP at two hidden layers and only lets `num_blocks` reach the residual
  net; here the one argument sets both.

- embedding_net:

  Optional summary network built with
  [`embedding_mlp()`](https://neuralsbi.pedrodelima.com/reference/embedding_mlp.md).
  The classifier then sees \\(\theta, f\_\psi(x))\\, with the embedding
  trained jointly. Ignored (with a warning) by `"logistic"`.

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

An object of class `nsbi_nre`. Evaluate the learned ratio with
[`log_ratio()`](https://neuralsbi.pedrodelima.com/reference/log_ratio.md)
or turn it into a posterior with
[`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md).

## Why a classifier instead of a density

A density estimator has to spend capacity describing the shape of a
distribution – normalization, tails, the lot – even where that shape
does not affect the posterior. A classifier only has to say which of two
parameter values explains the data better, and the optimal classifier's
logit *is* the log ratio. That makes NRE the natural choice when the
data are high-dimensional or awkward to model directly (discrete counts,
mixed types, anything a flow handles badly) but easy to discriminate.

Like [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) and
unlike [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md),
the ratio is learned for a single observation, so the log-likelihood of
\\n\\ independent trials from the same parameter is a sum:

\$\$\log \frac{p(x_1, \ldots, x_n \mid \theta)}{\prod_i p(x_i)} =
\sum\_{i=1}^{n} \log r(\theta, x_i).\$\$

Train once, condition on as many trials as you like. The price is the
same one [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md)
pays: posterior draws cost an MCMC run rather than a forward pass.

## The training objective

The default is the *atomic* loss of Durkan et al. (2020), which is what
`sbi`'s `NRE` (an alias for `NRE_B`) trains. For each simulation
\\(\theta_i, x_i)\\ in a minibatch, `num_atoms - 1` contrasting
parameters are taken from the other simulations in that batch, and the
classifier is scored on a `num_atoms`-way softmax over which of them
produced \\x_i\\:

\$\$\mathcal{L} = -\frac{1}{b}\sum\_{i} \left\[ f(\theta_i, x_i) - \log
\sum\_{k} \exp f(\theta\_{ik}, x_i) \right\].\$\$

More atoms mean a harder discrimination problem and a sharper ratio, at
a linear cost in forward passes per epoch. `num_atoms = 10` is `sbi`'s
default and this one. The softmax is invariant to adding any function of
\\x\\ to \\f\\, so the learned ratio is calibrated up to a constant at
each fixed observation – exactly what a posterior needs, and the reason
[`log_prob()`](https://neuralsbi.pedrodelima.com/reference/log_prob.md)
on an NRE posterior is unnormalized.

## Standardization has no Jacobian here

The estimators in
[`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) and
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) train on
z-scored data and need a change-of-variables term to report densities in
the original units. A ratio needs none: both \\p(x \mid \theta)\\ and
\\p(x)\\ pick up the same Jacobian factor and it cancels, so
[`log_ratio()`](https://neuralsbi.pedrodelima.com/reference/log_ratio.md)
in standardized space is already the ratio in the units the simulator
returned.

## References

Hermans, J., Begy, V. and Louppe, G. (2020). Likelihood-free MCMC with
Amortized Approximate Ratio Estimators. *ICML*.
[doi:10.48550/arXiv.1903.04057](https://doi.org/10.48550/arXiv.1903.04057)

Durkan, C., Murray, I. and Papamakarios, G. (2020). On Contrastive
Learning for Likelihood-free Inference. *ICML*.
[doi:10.48550/arXiv.2002.03712](https://doi.org/10.48550/arXiv.2002.03712)

## See also

[`log_ratio()`](https://neuralsbi.pedrodelima.com/reference/log_ratio.md)
to evaluate the ratio,
[`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md)
to sample it,
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) for the
likelihood factorization and
[`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) for the
posterior one.

## Examples

``` r
# One noisy measurement per simulator call; the observation is 50 of them.
prior <- prior_uniform(c(mu = -3), c(mu = 3))
simulator <- function(mu) c(y = rnorm(1, mu, 0.5))

fit <- nre(prior, simulator, n_simulations = 2000, classifier = "logistic")

x_obs <- matrix(rnorm(50, mean = 1, sd = 0.5), ncol = 1)
log_ratio(fit, theta = c(1), x = x_obs)
#> [1] 65.05644
```
