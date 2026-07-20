# Assemble batched lower-triangular Cholesky factors from the flat head output. Diagonal entries are passed through softplus (+ eps) to stay positive. Returns a tensor of shape (batch, K, p, p).

Assemble batched lower-triangular Cholesky factors from the flat head
output. Diagonal entries are passed through softplus (+ eps) to stay
positive. Returns a tensor of shape (batch, K, p, p).

## Usage

``` r
mdn_build_tril(net, tril_flat)
```
