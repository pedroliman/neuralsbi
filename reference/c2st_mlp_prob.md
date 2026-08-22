# `sbibm`'s C2ST classifier, trained on one fold

Adam on the binary log loss under the `MLPClassifier` settings `sbibm`
uses: minibatches of 200 in a fresh random order each epoch, learning
rate `1e-3`, and a stop once the epoch loss has failed to improve on the
best by more than `tol` for `n_iter_no_change` epochs in a row.

## Usage

``` r
c2st_mlp_prob(
  x_train,
  y_train,
  x_test,
  hidden,
  alpha = 1e-04,
  lr = 0.001,
  batch_size = 200L,
  max_epochs = 10000L,
  tol = 1e-04,
  n_iter_no_change = 10L,
  device = "cpu"
)
```

## Arguments

- x_train, y_train:

  Training draws and their 0/1 labels.

- x_test:

  Draws to score.

- hidden:

  Widths of the hidden layers, one entry per layer.

- alpha:

  L2 penalty on the weights.

- lr:

  Adam step size.

- batch_size:

  Minibatch size, capped at the number of training draws.

- max_epochs:

  Cap on epochs.

- tol:

  Smallest loss improvement that counts as progress.

- n_iter_no_change:

  Epochs without progress before training stops.

- device:

  Resolved torch device to train on.

## Value

Predicted probability of class 1 for each row of `x_test`.

## Details

`scikit-learn`'s `alpha` penalizes the weights and leaves the biases
alone, adding `alpha * W / batch` to the gradient. That is what
`weight_decay` does, so the optimizer carries it in two parameter groups
rather than the loss carrying it as a graph node: the penalty is also
part of the loss the stopping rule reads, but evaluating it once an
epoch instead of once a batch costs nothing and saves half the running
time of the loop.
