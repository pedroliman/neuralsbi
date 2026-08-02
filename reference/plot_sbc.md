# Plot an SBC rank histogram

Uniform bars indicate calibration; a U shape means the posterior is too
narrow (overconfident); an inverted-U means it is too wide.

## Usage

``` r
plot_sbc(sbc_result, param = 1L, bins = 20L)
```

## Arguments

- sbc_result:

  An `nsbi_sbc` object from
  [`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md).

- param:

  Which parameter to plot: an index (default 1), or a parameter name
  when [`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md) was
  run against a fit with named parameters (see
  [`prior_uniform()`](https://neuralsbi.pedrodelima.com/reference/prior_uniform.md)/[`prior_normal()`](https://neuralsbi.pedrodelima.com/reference/prior_normal.md)).
  The title uses that parameter's name when there is one, rendered as a
  plotmath symbol if the name parses as R syntax.

- bins:

  Number of histogram bins.

## Value

A `ggplot` object (also drawn as a side effect), invisibly.
