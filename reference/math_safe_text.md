# Vectorized, parse-safe plotmath label text

Re-quotes any entry that is not valid R syntax as a string literal,
which always parses (and renders as plain text under plotmath). Applying
this before a label is used as a facet/column name guarantees
[`ggplot2::label_parsed()`](https://ggplot2.tidyverse.org/reference/labellers.html)
never hits a parse error, whether or not the original label happens to
look like math (`"beta[1]"`) or not (`"growth rate"`).

## Usage

``` r
math_safe_text(labels)
```
