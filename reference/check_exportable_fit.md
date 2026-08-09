# Only an [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) fit has something to export

The generated Stan code is a transpiled *density*:
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md)'s
\\q\_\phi(x \mid \theta)\\ written out term by term. An
[`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md) fit holds
a classifier, not a density, and nothing in `R/stan.R` knows how to
write one out, so the two ways a user reaches this code say so by name
rather than failing on a bare
[`stopifnot()`](https://rdrr.io/r/base/stopifnot.html).

## Usage

``` r
check_exportable_fit(fit)
```
