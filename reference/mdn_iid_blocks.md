# Walk the observation chunks of an MDN's i.i.d. density, calling `collect(idx, lp)` with the `(n_theta, n_chunk)` tensor of log-densities

The MLP maps `theta` to the mixture parameters and never sees `x`, so
the forward pass and the Cholesky assembly run once for the whole
observation set and only the quadratic form is chunked.

## Usage

``` r
mdn_iid_blocks(de, xt, theta, max_batch, collect)
```
