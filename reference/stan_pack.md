# Flatten the trained weights into one vector, recording where each block sits

Stan has no closures and no struct types, so everything the generated
function needs arrives as a single `vector w` and is sliced back out
inside. Matrices are flattened column-major, which is exactly what
Stan's `to_matrix()` expects, so the two sides agree without any
transposition.

## Usage

``` r
stan_pack(fit)
```
