# Plot nominal vs. empirical credible-interval coverage

Well-calibrated posteriors lie on the diagonal. Curves above the
diagonal mean the posterior is too wide (conservative); below means
overconfident. A shaded band shows the Monte-Carlo uncertainty from the
finite number of SBC trials.

## Usage

``` r
plot_coverage(sbc_result, levels = seq(0.05, 0.95, by = 0.05))
```

## Arguments

- sbc_result:

  An `nsbi_sbc` object from
  [`sbc()`](https://pedroliman.github.io/neuralsbi/reference/sbc.md).

- levels:

  Nominal credibility levels to evaluate.

## Value

Invisibly, the coverage data frame from
[`expected_coverage()`](https://pedroliman.github.io/neuralsbi/reference/expected_coverage.md).
