# Dispatch to the requested ratio estimator

The counterpart of `fit_density_estimator()`. Arguments are forwarded
explicitly rather than through `...` and
[`formals()`](https://rdrr.io/r/base/formals.html), because there are
only two targets and they take disjoint arguments: the neural
classifiers take the whole training block, the closed-form one takes
none of it.

## Usage

``` r
fit_ratio_estimator(
  classifier,
  theta_z,
  x_z,
  hidden,
  n_blocks,
  num_atoms,
  embedding,
  max_epochs,
  batch_size,
  lr,
  validation_fraction,
  patience,
  n_restarts,
  clip_grad_norm,
  seed,
  verbose,
  device
)
```
