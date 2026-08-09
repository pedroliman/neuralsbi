# Run `expr` with the RNG parked at `seed`, restoring the stream afterwards

The caller's random stream is untouched, so a training run is still
reproducible from the seed it was given.

## Usage

``` r
with_fixed_seed(seed, expr)
```
