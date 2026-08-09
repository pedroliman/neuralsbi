# Quadratic feature basis for the logistic ratio estimator

`[1, z, vech(z z')]` for `z = (theta, x)`. The basis is chosen so the
estimator is *exact* for a linear-Gaussian simulator: there \\\log p(x
\mid \theta)\\ is a quadratic form in \\(\theta, x)\\, so the log
ratio's parameter dependence lies inside this span and the fit is
limited only by estimation error (see
[`fit_logistic_ratio()`](https://neuralsbi.pedrodelima.com/reference/fit_logistic_ratio.md)
for why the evidence term does not spoil that). It is the regression
oracle for [`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md)
that `"linear_gaussian"` is for
[`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) and
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md).

## Usage

``` r
nre_features(theta, x)
```

## Details

It costs `1 + d + d(d+1)/2` columns for `d = dim_theta + dim_x`, so it
is a baseline for small models, not a substitute for a neural classifier
on wide data.
