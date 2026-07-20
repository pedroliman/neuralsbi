# Expected coverage of central credible intervals

Uses the SBC ranks to compare nominal credible levels with the empirical
fraction of trials in which the true parameter falls inside the
corresponding central interval. Well-calibrated posteriors lie on the
diagonal.

## Usage

``` r
expected_coverage(sbc_result, levels = seq(0.05, 0.95, by = 0.05))
```

## Arguments

- sbc_result:

  An `nsbi_sbc` object from
  [`sbc()`](https://pedroliman.github.io/neural.sbi/reference/sbc.md).

- levels:

  Nominal credibility levels to evaluate.

## Value

A data frame with `nominal` and per-parameter empirical coverage.
