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
  verbose,
  device = "cpu"
)
```

## Details

`arch` carries the architecture fields specific to the caller –
`n_components`/`hidden` for the MDN, `n_transforms`/`hidden` for the
MAF, `n_transforms`/`hidden`/`n_bins`/`tail_bound` for the NSF – and is
spliced into the returned list ahead of `embedding`, matching the field
order each estimator returned before this helper existed. This helper
never needs to know what `arch`'s fields are.

`device` is a raw keyword here (`"cpu"`, `"cuda"`, `"mps"`, `"gpu"` or
`"auto"`) – resolving it to an actual, available device needs `torch`
loaded, so that happens inside
[`train_conditional_de()`](https://neuralsbi.pedrodelima.com/reference/train_conditional_de.md),
after its own argument checks
([`check_train_controls()`](https://neuralsbi.pedrodelima.com/reference/check_train_controls.md))
have already run without needing `torch` at all. The *resolved* string
comes back on
[`train_conditional_de()`](https://neuralsbi.pedrodelima.com/reference/train_conditional_de.md)'s
return value and is stored on the returned estimator (never a torch
device object, which would not survive
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html)) so
[`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md)/[`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md)
can see what it was actually fit with.
