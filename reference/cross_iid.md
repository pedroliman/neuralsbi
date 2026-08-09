# Walk the theta blocks of the cross product, calling `collect(idx, lp)`

`lp` is the flat vector `score()` returns for the block, running
theta-major: the whole observation set for the first parameter, then for
the second, and so on. `score` is
[de_log_prob()](https://neuralsbi.pedrodelima.com/reference/density_estimator.md)
for a density estimator and
[`nre_score()`](https://neuralsbi.pedrodelima.com/reference/nre_score.md)
for a ratio estimator; the blocking is the same either way.

## Usage

``` r
cross_iid(de, x, theta, max_batch, collect, score = de_log_prob)
```
