# Independent RNG streams, one per chunk

Derived from the caller's current RNG state, so they are reproducible
under [`set.seed()`](https://rdrr.io/r/base/Random.html) yet independent
of the backend. Consumes exactly one draw from the caller's stream and
leaves its kind untouched.

## Usage

``` r
rng_streams(n)
```
