# Train a neural ratio estimator on standardized (theta, x)

Passes `min_val_rows = 2L` down to
[`fit_torch_de()`](https://neuralsbi.pedrodelima.com/reference/fit_torch_de.md),
unlike the MDN/MAF/NSF callers. The atomic objective
([`nre_atomic_log_prob()`](https://neuralsbi.pedrodelima.com/reference/nre_atomic_log_prob.md))
needs a second row to contrast the true parameter against; with a
validation split of one row it silently returns a constant zero loss
every epoch instead of a real signal, which breaks early stopping
(GitHub \#188).

## Usage

``` r
fit_nre_net(
  theta,
  x,
  classifier = "resnet",
  hidden = 50L,
  n_blocks = 2L,
  num_atoms = 10L,
  max_epochs = 2000L,
  batch_size = 200L,
  lr = 5e-04,
  validation_fraction = 0.1,
  patience = 20L,
  n_restarts = 1L,
  clip_grad_norm = 5,
  embedding = NULL,
  seed = NULL,
  verbose = FALSE,
  device = "cpu"
)
```
