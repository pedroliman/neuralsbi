# Rebuild an estimator's network from the architecture recorded on the fit

The one place that knows how to turn a stored estimator back into a
torch module. Every field it reads is set by the matching `fit_*()`, so
adding an estimator means adding a branch here.

## Usage

``` r
de_rebuild_net(de)
```
