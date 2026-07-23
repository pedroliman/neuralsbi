# Apply one spline transform elementwise over all theta dimensions.

Spline parameters are computed from `z` (and `x`); the transform itself
is applied to `values` (defaults to `z`). The separation matters when
inverting: parameters must come from the partially reconstructed theta
while the inverse acts on the base-space values.

## Usage

``` r
nsf_apply(net, made, z, x, inverse = FALSE, values = z)
```
