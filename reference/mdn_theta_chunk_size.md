# How many rows of `theta` to run through the MLP/Cholesky assembly at once

Mirrors the theta-blocking
[`cross_iid()`](https://neuralsbi.pedrodelima.com/reference/cross_iid.md)
does for the flow estimators, so an MDN bounds the `(theta, x)` pair
count the same way:
[`mdn_mixture()`](https://neuralsbi.pedrodelima.com/reference/mdn_mixture.md)
materializes a `(theta_chunk, K, dim_theta, dim_theta)` tensor for one
block, and `theta_chunk * n_obs <= max_batch` keeps that bounded
regardless of which of `theta` or `x` is the large dimension (#240).

## Usage

``` r
mdn_theta_chunk_size(n_theta, n_obs, max_batch)
```
