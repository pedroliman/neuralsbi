# Stack of `(linear, relu)` pairs, masked or plain.

Each consecutive pair in `dims` becomes one `(linear, relu)` step,
matching the trunk-building loop shared by the MADE, MDN and embedding
networks.

## Usage

``` r
mlp_layers(dims, masks = NULL)
```

## Arguments

- dims:

  Layer sizes in order, e.g. `c(input, hidden1, hidden2, ...)`.

- masks:

  One weight mask per pair (as returned by `made_masks()$hidden`) for an
  autoregressive trunk; `NULL` for a plain MLP.

## Value

A list of torch layers, ready for `do.call(torch::nn_sequential, .)`.
