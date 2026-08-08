# Shared training engine for torch conditional density estimators

All neural estimators (MDN, MAF, NSF) share one training loop so that
robustness features are implemented once: train/validation split, Adam,
minibatching, early stopping on validation loss, learning-rate decay on
plateau, gradient clipping, and best-of-`n_restarts` reinitialization.
The defaults (batch 200, lr 5e-4, 10% validation, patience 20, clip norm
5) match Python `sbi`, so results are comparable across the two
packages.

## Usage

``` r
train_conditional_de(
  build_net,
  log_prob_fn,
  theta,
  x,
  max_epochs = 2000L,
  batch_size = 200L,
  lr = 5e-04,
  validation_fraction = 0.1,
  patience = 20L,
  n_restarts = 1L,
  clip_grad_norm = 5,
  lr_patience = 10L,
  lr_factor = 0.5,
  min_lr = 1e-06,
  seed = NULL,
  verbose = FALSE,
  device = "cpu"
)
```

## Arguments

- build_net:

  A zero-argument function returning a *fresh* torch module. Called once
  per restart so each restart gets new initial weights.

- log_prob_fn:

  `function(net, theta, x)` returning a length-`b` tensor of log
  densities for a `(b, p)` theta tensor and `(b, q)` x tensor.

- theta, x:

  Standardized training matrices.

- n_restarts:

  Train this many independently initialized networks and keep the one
  with the best validation loss.

- clip_grad_norm:

  Maximum gradient norm (set `Inf` to disable).

- lr_patience, lr_factor, min_lr:

  Reduce the learning rate by `lr_factor` after `lr_patience` epochs
  without validation improvement, down to `min_lr`.

- device:

  Torch device keyword: `"cpu"` (the default), `"cuda"`, `"mps"`, or
  `"gpu"`/`"auto"` (see `resolve_device()`, which is what turns this
  into an actual, available device). Training and validation tensors are
  created there, and the net is moved there right after `build_net()`,
  so the two never disagree the way they do under a bare
  [`torch::with_device()`](https://torch.mlverse.org/docs/reference/local_device.html).

## Value

`list(net, best_val_loss, history, device)`, where `history` is a data
frame of per-epoch train/validation losses for the winning restart and
`device` is the resolved device (`"cpu"`, `"cuda"` or `"mps"`) training
actually ran on.
