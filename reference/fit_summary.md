# Print a fit and return its training metadata invisibly

The body every [`summary()`](https://rdrr.io/r/base/summary.html) method
for a fit shares. `estimator` is the one field whose name depends on
what was trained: `density_estimator` for
[`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) and
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md),
`classifier` for
[`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md).
[`print()`](https://rdrr.io/r/base/print.html) inside dispatches on the
object's own class, so nothing else here has to know which fit it has.

## Usage

``` r
fit_summary(object, estimator)
```

## Arguments

- object:

  The fit.

- estimator:

  A one-element named list, spliced in ahead of the rest.
