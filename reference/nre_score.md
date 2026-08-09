# The scorer [`cross_iid()`](https://neuralsbi.pedrodelima.com/reference/cross_iid.md) needs, with the ratio's argument order

[`cross_iid()`](https://neuralsbi.pedrodelima.com/reference/cross_iid.md)
calls `score(de, target, condition)` – data first, parameters second,
the order
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md)'s
estimator is trained in. A ratio has no target and no condition, and
reads better with the parameter first, so the flip happens here rather
than in
[`de_log_ratio()`](https://neuralsbi.pedrodelima.com/reference/de_log_ratio.md)'s
signature.

## Usage

``` r
nre_score(de, x, theta)
```
