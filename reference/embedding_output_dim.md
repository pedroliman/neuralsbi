# Effective conditioning dimension after an (optional) embedding.

The identity embedding (`spec = NULL`) leaves the data dimension
unchanged; otherwise the estimator conditions on `output_dim` features.

## Usage

``` r
embedding_output_dim(spec, dim_x)
```
