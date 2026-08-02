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
  [`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md).

- levels:

  Nominal credibility levels to evaluate, each strictly between 0 and 1.
  Passed to
  [`expected_coverage()`](https://neuralsbi.pedrodelima.com/reference/expected_coverage.md),
  which checks them.

## Value

A `ggplot` object (also drawn as a side effect), invisibly.
