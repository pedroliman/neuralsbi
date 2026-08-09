# The atomic NRE-B objective, in the shape [`train_conditional_de()`](https://neuralsbi.pedrodelima.com/reference/train_conditional_de.md) wants

The training engine minimizes `-mean(log_prob_fn(net, theta, x))`, so
the per-row quantity returned here is the atomic log-softmax: the true
parameter's logit minus the log-sum-exp over its atoms. Maximizing it is
minimizing the cross-entropy of picking the true parameter out of the
`num_atoms` candidates.

## Usage

``` r
nre_atomic_log_prob(num_atoms)
```

## Details

`num_atoms` is clamped to the batch size, as `sbi` does, so a short
final minibatch or a small validation split does not ask for more
contrasts than the batch can supply. A batch of one has no contrast at
all and contributes nothing, but it still has to contribute a tensor
torch built: a constant of its own making has no graph for `backward()`
to walk.
[`minibatches()`](https://neuralsbi.pedrodelima.com/reference/minibatches.md)
keeps single-row batches out of training; a validation split of one row
still arrives here.
