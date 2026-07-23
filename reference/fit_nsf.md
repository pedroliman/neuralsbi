# Train an NSF on standardized (theta, x)

Train an NSF on standardized (theta, x)

## Usage

``` r
fit_nsf(
  theta,
  x,
  n_transforms = 5L,
  hidden = c(50L, 50L),
  n_bins = 10L,
  tail_bound = 3,
  max_epochs = 2000L,
  batch_size = 200L,
  lr = 5e-04,
  validation_fraction = 0.1,
  patience = 20L,
  n_restarts = 1L,
  clip_grad_norm = 5,
  embedding = NULL,
  seed = NULL,
  verbose = FALSE
)
```
