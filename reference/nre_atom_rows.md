# Which rows of a minibatch each simulation is scored against

Returns a length-`b * k` vector of row indices, `k` per simulation and
running simulation-major: the true parameter for row `i` first, then
`k - 1` contrasts drawn without replacement from the *other* rows of the
batch. That is `sbi`'s atom construction, done with R's RNG rather than
`torch.multinomial` so it is seeded by the same
[`set.seed()`](https://rdrr.io/r/base/Random.html) the training loop's
split and batch order already are.

## Usage

``` r
nre_atom_rows(b, k, deterministic = FALSE)
```

## Details

`deterministic` freezes the draw. The atoms are resampled every time the
loss is evaluated, which is the objective's business during training but
makes a poor early-stopping signal: the validation loss would move
between epochs because the contrasts changed, not because the classifier
did. Freezing it for the validation pass costs nothing and makes the two
numbers comparable. (`sbi` resamples there too, and pays for it in
noisier stopping.)
