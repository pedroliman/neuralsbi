# Validate the training controls

Every one of these decides how the split, the batches or the restart
loop is built, and an unchecked value is reported by whichever base
function hits it first: `batch_size = 0` comes back as "invalid '(to -
from)/by'", `validation_fraction = 1` as "wrong sign in 'by' argument",
and `n_restarts = 0` as "Training failed: no restart produced a finite
validation loss", which blames training for an argument.
[`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) and
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) call this
before they simulate, so a typo does not cost the budget first.

## Usage

``` r
check_train_controls(
  max_epochs,
  batch_size,
  lr,
  validation_fraction,
  patience,
  n_restarts,
  clip_grad_norm,
  n = NULL
)
```

## Arguments

- max_epochs, batch_size, lr, validation_fraction, patience:

  Neural training controls (Adam optimizer, early stopping on validation
  loss). The defaults (`batch_size = 200`, `lr = 5e-4`,
  `validation_fraction = 0.1`, `patience = 20`) match Python `sbi`;
  `max_epochs` is a high guard cap that early stopping normally reaches
  first.

- n_restarts:

  Train this many independently initialized networks and keep the one
  with the best validation loss (guards against bad initializations and
  MDN mode collapse).

- clip_grad_norm:

  Maximum gradient norm during training (`Inf` disables clipping). The
  learning rate also decays 2x after 10 epochs without validation
  improvement.

- n:

  Number of training rows, or `NULL` when they do not exist yet.

## Details

`n` is optional because that call happens before there are any rows.
When it is known, `validation_fraction` is checked against it: the
requirement is that both sides of the split come out non-empty, which
the fraction alone cannot decide.
