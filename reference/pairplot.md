# Visualize posterior samples

A dependency-free (base graphics) pair plot: 1-D marginal densities on
the diagonal and 2-D scatter/contours off-diagonal, with optional
markers for a reference (e.g. true) parameter value. Analogous to
`sbi`'s `pairplot`.

## Usage

``` r
pairplot(
  samples,
  truth = NULL,
  labels = NULL,
  limits = NULL,
  col = grDevices::adjustcolor("steelblue", 0.4),
  ...
)
```

## Arguments

- samples:

  A matrix of posterior draws (rows = draws), or an `nsbi_samples`
  object.

- truth:

  Optional reference parameter vector to overlay.

- labels:

  Optional parameter labels.

- limits:

  Optional list/matrix of per-parameter c(lo, hi) axis limits.

- col:

  Point colour.

- ...:

  Passed to plotting calls.

## Value

Invisibly, the samples.
