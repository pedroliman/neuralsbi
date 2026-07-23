# Apply an estimator's embedding to conditioning data, if it has one.

Estimators call this once per forward/inverse pass so the embedding runs
a single time; the raw (standardized) `x` still enters at the `de_*`
boundary, keeping `de$dim_x` the raw data dimension.

## Usage

``` r
embed_x(net, x)
```
