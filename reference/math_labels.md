# Vectorized, parse-safe plotmath labels

Like
[`math_expr()`](https://neuralsbi.pedrodelima.com/reference/math_expr.md)
but for a whole vector at once, returning an
[`expression()`](https://rdrr.io/r/base/expression.html) (so a discrete
ggplot2 scale's `labels =` can render each entry as plotmath, mirroring
[`ggplot2::label_parsed()`](https://ggplot2.tidyverse.org/reference/labellers.html)
for facet strips).

## Usage

``` r
math_labels(labels)
```
