# Train an MDN on standardized (theta, x)

Train an MDN on standardized (theta, x)

## Usage

``` r
fit_mdn(
  theta,
  x,
  n_components = 10L,
  hidden = c(50L, 50L),
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

## Arguments

- embedding:

  Optional embedding-network spec (see
  [`embedding_mlp()`](https://neuralsbi.pedrodelima.com/reference/embedding_mlp.md));
  the MDN then conditions on the learned features instead of the raw
  `x`.
