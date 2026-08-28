# Walk the theta and observation blocks of an MDN's i.i.d. density, calling `collect(theta_idx, obs_idx, lp)` with the `(n_theta_chunk, n_obs_chunk)` tensor of log-densities

The MLP maps `theta` to the mixture parameters and never sees `x`, so
within one theta block the forward pass and the Cholesky assembly run
once and only the quadratic form is chunked over observations. `theta`
itself is blocked first
([`mdn_theta_chunk_size()`](https://neuralsbi.pedrodelima.com/reference/mdn_theta_chunk_size.md)),
so `max_batch` bounds
[`mdn_mixture()`](https://neuralsbi.pedrodelima.com/reference/mdn_mixture.md)'s
output the same way
[`cross_iid()`](https://neuralsbi.pedrodelima.com/reference/cross_iid.md)
bounds a flow's forward pass, rather than only chunking the observation
side (#240).

## Usage

``` r
mdn_iid_blocks(de, xt, theta, max_batch, collect)
```
