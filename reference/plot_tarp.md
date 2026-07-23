# Plot TARP expected coverage

Draws the expected coverage probability (ECP) curve from
[`tarp()`](https://pedroliman.github.io/neuralsbi/reference/tarp.md)
against the nominal credibility level. A calibrated posterior lies on
the diagonal; a curve above the diagonal means the posterior is too wide
(conservative), below means overconfident. The shaded band shows the
Monte-Carlo uncertainty from the finite number of TARP trials.

## Usage

``` r
plot_tarp(tarp_result)
```

## Arguments

- tarp_result:

  An `nsbi_tarp` object from
  [`tarp()`](https://pedroliman.github.io/neuralsbi/reference/tarp.md).

## Value

Invisibly, a data frame with `nominal` and `ecp` columns.
