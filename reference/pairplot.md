# Visualize posterior samples

A pair plot built on
[`GGally::ggpairs()`](https://ggobi.github.io/ggally/reference/ggpairs.html):
1-D marginal densities on the diagonal and 2-D highest-density regions
(via
[`ggdensity::geom_hdr()`](https://jamesotto852.github.io/ggdensity/reference/geom_hdr.html))
in the lower triangle, with optional markers for a reference (e.g. true)
parameter value. Analogous to Python `sbi`'s `pairplot`.

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
  matrix of per-parameter axis limits. Defaults to each parameter's own
  data range (padded 5%), applied consistently across every panel that
  plots it – without this, the lower-triangle and diagonal panels for
  the same parameter can each draw a different range, since
  [`ggdensity::geom_hdr()`](https://jamesotto852.github.io/ggdensity/reference/geom_hdr.html)
  and
  [`ggplot2::geom_density()`](https://ggplot2.tidyverse.org/reference/geom_density.html)
  estimate their density grid independently per panel.

- col:

  Density-region and marginal-density fill colour.

- alpha:

  Marginal-density fill transparency. The lower-triangle highest-density
  regions shade themselves by probability level (99/95/80/50%) instead,
  via `ggdensity`'s own `alpha` mapping.

- ...:

  Passed to the lower-triangle
  [`ggdensity::geom_hdr()`](https://jamesotto852.github.io/ggdensity/reference/geom_hdr.html)
  layer.

## Value

A `ggmatrix` object (also drawn as a side effect), invisibly.
