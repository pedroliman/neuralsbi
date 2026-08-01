# Parse a label as a plotmath expression when possible

Parameter/outcome names that happen to be valid R syntax (`"beta[1]"`,
`"rho"`, `"sigma^2"`) render as their mathematical symbol – Greek
letters, sub/superscripts – when passed through R's plotmath. Names that
are not parseable, or `NULL`, pass through unchanged so callers can fall
back to a plain-text label.

## Usage

``` r
math_expr(label)
```
