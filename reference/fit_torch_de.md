# Shared body behind `fit_mdn()`, `fit_maf()` and `fit_nsf()`

Coerces `theta` and `x`, builds the net from
`build_net_fn(dim_x, dim_theta)` now that both are known, trains it with
[`train_conditional_de()`](https://neuralsbi.pedrodelima.com/reference/train_conditional_de.md),
and packages the result into a fitted `nsbi_de` object.

## Usage

``` r
fit_torch_de(
  theta,
  x,
  build_net_fn,
  log_prob_fn,
  class,
  arch,
  max_epochs,
  batch_size,
  lr,
  validation_fraction,
  patience,
  n_restarts,
  clip_grad_norm,
  embedding,
  seed,
  verbose
)
```

## Details

`arch` carries the architecture fields specific to the caller –
`n_components`/`hidden` for the MDN, `n_transforms`/`hidden` for the
MAF, `n_transforms`/`hidden`/`n_bins`/`tail_bound` for the NSF – and is
spliced into the returned list ahead of `embedding`, matching the field
order each estimator returned before this helper existed. This helper
never needs to know what `arch`'s fields are.
