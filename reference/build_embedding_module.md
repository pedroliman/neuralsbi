# Build the embedding torch submodule, or `NULL` for the identity embedding.

Constructed lazily so no torch object exists at package-load time.
Returns an instantiated `nn_module` (call site stores it as a submodule
so its parameters train jointly and travel with the estimator's
`state_dict`).

## Usage

``` r
build_embedding_module(spec, dim_x)
```
