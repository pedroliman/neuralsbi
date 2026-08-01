# Walk the theta blocks of the cross product, calling `collect(idx, lp)`

`lp` is the flat vector
[de_log_prob()](https://neuralsbi.pedrodelima.com/reference/density_estimator.md)
returns for the block, running theta-major: the whole observation set
for the first parameter, then for the second, and so on.

## Usage

``` r
cross_iid(de, x, theta, max_batch, collect)
```
