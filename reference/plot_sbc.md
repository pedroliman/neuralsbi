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
  [`sbc()`](https://pedroliman.github.io/neuralsbi/reference/sbc.md).

- param:

  Which parameter index to plot (default 1). The title uses that
  parameter's name (from a named prior, see
  [`prior_uniform()`](https://pedroliman.github.io/neuralsbi/reference/prior_uniform.md)/
  [`prior_normal()`](https://pedroliman.github.io/neuralsbi/reference/prior_normal.md))
  when
  [`sbc()`](https://pedroliman.github.io/neuralsbi/reference/sbc.md) was
  run against a fit that has one, rendered as a plotmath symbol if the
  name parses as R syntax.

- bins:

  Number of histogram bins.

## Value

A `ggplot` object (also drawn as a side effect), invisibly.
