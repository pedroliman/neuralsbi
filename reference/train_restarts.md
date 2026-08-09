# Restart loop behind `train_conditional_de()`

Split out so
[`train_conditional_de()`](https://neuralsbi.pedrodelima.com/reference/train_conditional_de.md)
can wrap the whole loop – all restarts, one progress bar – in a
progress-reporting context.

## Usage

``` r
train_restarts(
  build_net,
  log_prob_fn,
  theta,
  x,
  max_epochs,
  batch_size,
  lr,
  validation_fraction,
  patience,
  n_restarts,
  clip_grad_norm,
  lr_patience,
  lr_factor,
  min_lr,
  seed,
  verbose,
  device = "cpu"
)
```

## Details

`device` here is still the raw, unresolved keyword
[`train_conditional_de()`](https://neuralsbi.pedrodelima.com/reference/train_conditional_de.md)
was given (`"cpu"`, `"cuda"`, `"mps"`, `"gpu"` or `"auto"`); resolving
it to an actual, available device needs `torch` loaded (see
[`resolve_device()`](https://neuralsbi.pedrodelima.com/reference/resolve_device.md)),
and this is the first point that is guaranteed true –
[`require_torch()`](https://neuralsbi.pedrodelima.com/reference/require_torch.md)
is the line above. Doing it here rather than earlier in
[`train_conditional_de()`](https://neuralsbi.pedrodelima.com/reference/train_conditional_de.md)
keeps
[`check_train_controls()`](https://neuralsbi.pedrodelima.com/reference/check_train_controls.md)
(which needs no torch at all) running first, so a bad `batch_size` is
still reported before an unavailable device is.
