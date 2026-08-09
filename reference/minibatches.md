# Split a shuffled row order into minibatches

Stepping through the rows by `batch_size` leaves a short final batch,
and a final batch of *one* row is worse than short.
[`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md)'s atomic
loss has no contrast to score a single simulation against, so it returns
a constant that torch did not build and `backward()` errors with
"element 0 of tensors does not require grad"; every other estimator
takes a gradient step from one sample. Those rows join the previous
batch instead. It costs one batch of `batch_size + 1` once an epoch and
removes a failure that depends on nothing but `n_simulations` modulo
`batch_size`.

## Usage

``` r
minibatches(order, batch_size)
```

## Arguments

- order:

  Row indices, already shuffled.

- batch_size:

  Rows per batch.

## Value

A list of index vectors covering `order` in order.
