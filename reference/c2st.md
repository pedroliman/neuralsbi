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

  Matrices of samples (rows = draws, cols = dimensions).

- n_folds:

  Number of cross-validation folds.

- seed:

  Optional seed.

## Value

A list with mean CV accuracy and per-fold accuracies.
