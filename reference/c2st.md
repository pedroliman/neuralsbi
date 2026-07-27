# Classifier two-sample test (C2ST)

Trains a logistic-regression classifier to distinguish samples in `x`
from samples in `y` using cross-validation. A test accuracy near 0.5
means the two sample sets are indistinguishable (good); near 1.0 means
they differ. This is the standard SBI metric for comparing an estimated
posterior to a reference (e.g. an analytic posterior or long-run MCMC
draws).

## Usage

``` r
c2st(x, y, n_folds = 5L, seed = NULL)
```

## Arguments

- x, y:

  Matrices of samples (rows = draws, cols = dimensions). Sizes need not
  match; the larger is subsampled down to the smaller.

- n_folds:

  Number of cross-validation folds.

- seed:

  Optional seed.

## Value

A list with mean CV accuracy and per-fold accuracies.

## Details

The classifier is linear, which decides what the number can and cannot
see. It picks up a shift in location readily; it is close to blind to
two sample sets that share a mean and differ in spread or in the shape
of their dependence, because no hyperplane separates those. Python `sbi`
uses an MLP, so its C2ST sees more, and a 0.5 from here is the weaker
claim of the two. Read it alongside the moments rather than on its own.

Unequal sample sizes are balanced by subsampling the larger set, because
accuracy against unbalanced classes is not a two-sample test: 8000 draws
against 2000 identical ones scores 0.8 for a classifier that has learned
nothing except to always answer with the bigger class.
