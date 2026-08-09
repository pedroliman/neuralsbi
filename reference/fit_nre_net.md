# Train a neural ratio estimator on standardized (theta, x)

Train a neural ratio estimator on standardized (theta, x)

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
