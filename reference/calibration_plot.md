# Shared skeleton for a nominal-vs-empirical calibration plot

The Monte-Carlo band, the diagonal reference line, the equal-aspect
`[0, 1] x [0, 1]` frame and the minimal theme are the same figure in
[`plot_coverage()`](https://neuralsbi.pedrodelima.com/reference/plot_coverage.md)
and
[`plot_tarp()`](https://neuralsbi.pedrodelima.com/reference/plot_tarp.md);
only the curve drawn on top (one per parameter, with a legend, for
[`plot_coverage()`](https://neuralsbi.pedrodelima.com/reference/plot_coverage.md);
a single curve for
[`plot_tarp()`](https://neuralsbi.pedrodelima.com/reference/plot_tarp.md))
and the title differ. This builds the shared part and leaves the caller
to add its own
[`ggplot2::geom_line()`](https://ggplot2.tidyverse.org/reference/geom_path.html)
and
[`ggplot2::labs()`](https://ggplot2.tidyverse.org/reference/labs.html)
on top.

## Usage

``` r
calibration_plot(df, band, xlab, ylab)
```

## Arguments

- df:

  Base data for the plot, inherited by whatever curve layer the caller
  adds afterward.

- band:

  A data frame with `nominal`, `lo` and `hi` columns (see
  [`binom_band()`](https://neuralsbi.pedrodelima.com/reference/binom_band.md)),
  shaded as the Monte-Carlo uncertainty ribbon.

- xlab, ylab:

  Axis labels.

## Value

A `ggplot` object carrying the shared layers, ready for the caller's
curve, scale and title.
