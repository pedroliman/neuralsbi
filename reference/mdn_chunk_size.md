# How many observations to score in one call

`max_batch` counts (theta, x) pairs, so that is what sizes the chunk.
The quadratic form materializes `K` components and `p` dimensions per
pair on top of that, so a second cap keeps the `(T, K, n, p)`
intermediate to a few tens of megabytes however many components the
mixture has. Dividing by `K * p` alone, as this used to, chunked a
5000-observation call ten ways and paid the fixed cost of a torch call
ten times for a 400 KB intermediate.

## Usage

``` r
mdn_chunk_size(n_theta, max_batch, per_pair)
```
