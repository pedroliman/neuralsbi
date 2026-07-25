# Visualize posterior samples

A pair plot built on
[`GGally::ggpairs()`](https://ggobi.github.io/ggally/reference/ggpairs.html):
1-D marginal densities on the diagonal and 2-D scatter in the lower
triangle, with optional markers for a reference (e.g. true) parameter
value. Analogous to Python `sbi`'s `pairplot`.

## Usage

``` r
pairplot(
  samples,
  truth = NULL,
  labels = NULL,
  limits = NULL,
  col = "steelblue",
  alpha = 0.4,
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

  Optional parameter labels. Defaults to `colnames(samples)` (set
  automatically for
  [`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md)
  draws from a fit with named parameters – see
  [`prior_uniform()`](https://neuralsbi.pedrodelima.com/reference/prior_uniform.md)/[`prior_normal()`](https://neuralsbi.pedrodelima.com/reference/prior_normal.md))
  or `theta[1]`, `theta[2]`, .... Labels that parse as R syntax
  (`"beta[1]"`, `"rho"`) render as their plotmath symbol.

- limits:

  Optional list (one `c(lo, hi)` per parameter, in column order) or
  matrix of per-parameter axis limits.

- col:

  Point and density fill colour.

- alpha:

  Point and density fill transparency.

- ...:

  Passed to the lower-triangle
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)
  layer.

## Value

A `ggmatrix` object (also drawn as a side effect), invisibly.
