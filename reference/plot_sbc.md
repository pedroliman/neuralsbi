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
  [`sbc()`](https://pedroliman.github.io/neural.sbi/reference/sbc.md).

- param:

  Which parameter index to plot (default 1).

- bins:

  Number of histogram bins.

## Value

Invisibly, the rank vector.
