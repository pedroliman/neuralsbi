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
  verbose
)
```
