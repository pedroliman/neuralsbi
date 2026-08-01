# Report progress for one phase of work

Wraps a phase (simulation, training) so the built-in bar or a progressr
handler is active for its duration. `expr` is passed on unevaluated so
progressr can catch the progression conditions it signals.

## Usage

``` r
with_nsbi_progress(expr)
```
