# Evaluate `expr` under a given L'Ecuyer-CMRG stream

The sequential counterpart of `future(..., seed = stream)`, so
sequential and parallel runs draw the same random numbers.

## Usage

``` r
with_rng_stream(stream, expr)
```
