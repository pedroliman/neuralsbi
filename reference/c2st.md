# Classifier two-sample test (C2ST)

Trains a classifier to tell the draws in `x` apart from the draws in `y`
and reports its cross-validated test accuracy. An accuracy near 0.5
means the two sample sets are indistinguishable (good); near 1.0 means
they differ. This is the headline metric of the `sbibm` benchmark suite,
used there to score an approximate posterior against reference draws.

## Usage

``` r
c2st(
  x,
  y,
  n_folds = 5L,
  seed = NULL,
  classifier = c("mlp", "logistic"),
  z_score = TRUE,
  noise_scale = NULL,
  hidden = NULL,
  max_epochs = 10000L,
  device = "cpu"
)
```

## Arguments

- x, y:

  Matrices of samples (rows = draws, cols = dimensions). The two must be
  the same width, since they are draws of the same quantity. Row counts
  need not match; the larger set is subsampled down to the smaller. Pass
  reference draws as `x`, following `sbibm`.

- n_folds:

  Number of cross-validation folds. At least 2, and fewer than the
  number of draws in the smaller sample set.

- seed:

  Optional seed. Fixes the fold split, the subsampling, the network
  initialization and the minibatch order.

- classifier:

  `"mlp"` for the `sbibm` network, `"logistic"` for logistic regression.

- z_score:

  Standardize both sample sets by the mean and standard deviation of `x`
  before training. On by default, as in `sbibm`.

- noise_scale:

  Standard deviation of Gaussian noise added to both sample sets after
  z-scoring. `NULL`, the default, adds none. Set it when the draws are
  discrete or lie on a lower-dimensional set, where a classifier can
  separate the two sides on an artefact of representation.

- hidden:

  Widths of the hidden layers, one entry per layer. `NULL`, the default,
  is `sbibm`'s two layers of `10 * d` units.

- max_epochs:

  Cap on training epochs per fold. A guard, not a budget: the
  no-improvement rule normally stops training first.

- device:

  Torch device to train the network on: `"cpu"` (the default), `"cuda"`,
  `"mps"`, or `"gpu"`/`"auto"` for whichever is available. Ignored by
  `classifier = "logistic"`.

## Value

A list with the mean cross-validated accuracy and ROC AUC, the per-fold
values of each, the number of draws per class, the classifier used, and
a one-line reading of the accuracy.

## Details

The defaults reproduce the procedure in `sbibm/metrics/c2st.py`, so a
number from here is comparable with a published one: z-score both sample
sets using the mean and standard deviation of `x`, fit a
two-hidden-layer ReLU network of `10 * d` units per layer (`d` = number
of columns) with Adam, L2 penalty `1e-4`, minibatches of 200 and a
learning rate of `1e-3`, and average the accuracy over 5 shuffled
cross-validation folds. Training stops when the epoch loss fails to
improve by `1e-4` for 10 epochs in a row, which is what ends the fit
long before `max_epochs`. `sbibm` passes reference draws first, so pass
the reference as `x`: that is the set whose moments set the scale.

The network trains on torch, like every other neural piece of this
package, so `classifier = "mlp"` needs torch installed and errors
without it. That is the only part of the diagnostics that does.

Budget for the fit. `sbibm`'s own comparison, 10000 draws against 10000
in two dimensions, takes about 25 seconds on one CPU core; at `d = 5` it
is several minutes, because the hidden layers are `10 * d` wide and the
loss takes longer to flatten. Five folds of an all but unbounded epoch
budget is what the `sbibm` recipe asks for, and that is what it costs.
`hidden` and `max_epochs` cut it down, and `device` moves the fit to a
GPU.

`classifier = "logistic"` swaps the network for cross-validated logistic
regression. It needs no torch, it answers in milliseconds, and it is
what this function used before it was aligned with `sbibm`. It is also
linear: it sees a shift in location and is close to blind to two sample
sets that share a mean and differ in spread or in the shape of their
dependence. Use it as a cheap screen, or in a session with no torch, and
report the MLP number.

Unequal sample sizes are balanced by subsampling the larger set, which
`sbibm` does not do because it always compares 10000 draws against
10000. Accuracy against unbalanced classes is not a two-sample test:
8000 draws against 2000 identical ones scores 0.8 for a classifier that
has learned nothing except to always answer with the bigger class.

## References

Lopez-Paz, D. and Oquab, M. (2017). Revisiting classifier two-sample
tests. *ICLR*.
[doi:10.48550/arXiv.1610.06545](https://doi.org/10.48550/arXiv.1610.06545)

Lueckmann, J.-M., Boelts, J., Greenberg, D. S., Goncalves, P. J. and
Macke, J. H. (2021). Benchmarking simulation-based inference. *AISTATS*.
[doi:10.48550/arXiv.2101.04653](https://doi.org/10.48550/arXiv.2101.04653)

## Examples

``` r
a <- matrix(rnorm(400), ncol = 2)
b <- matrix(rnorm(400), ncol = 2)
c2st(a, b, classifier = "logistic", seed = 1)$accuracy
#> [1] 0.5175
```
