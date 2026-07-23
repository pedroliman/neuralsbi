# Embedding (summary) networks for structured observations

Raw observations are often high-dimensional or structured (a time
series, a set of summary statistics, an image) where feeding `x`
straight into the density estimator wastes capacity. An embedding
network learns a low- dimensional summary \\h = f\_\psi(x)\\ jointly
with the density estimator, so the conditioning path becomes
\\q\_\phi(\theta \mid f\_\psi(x))\\. This mirrors `sbi`'s
`embedding_net` argument.

## Usage

``` r
embedding_mlp(output_dim = 16L, hidden = c(64L, 64L))
```

## Arguments

- output_dim:

  Number of summary features the network emits. This is the effective
  data dimension the density estimator conditions on.

- hidden:

  Integer vector of hidden-layer widths (ReLU between layers). An empty
  vector gives a single linear map to `output_dim`.

## Value

An `nsbi_embedding` specification. It carries no torch objects (the
network is built lazily at fit time), so it is safe to construct without
`torch` installed.

## Details

`embedding_mlp()` builds a multilayer-perceptron summary network: a
stack of fully connected ReLU layers mapping the (standardized) data to
a vector of `output_dim` features. Pass the result to
[`npe()`](https://pedroliman.github.io/neural.sbi/reference/npe.md) via
`embedding_net`; it is trained end to end with the estimator and its
parameters live inside the fitted network, so sampling and `log_prob`
route through it automatically.

The embedding consumes the standardized data (the same z-scoring
[`npe()`](https://pedroliman.github.io/neural.sbi/reference/npe.md)
applies to `x` without an embedding), which keeps the summary network's
inputs on a common scale. Standardization of the *features* is
intentionally left to the network itself; the estimators operate on the
raw embedding output.

## See also

[`npe()`](https://pedroliman.github.io/neural.sbi/reference/npe.md)

## Examples

``` r
emb <- embedding_mlp(output_dim = 8, hidden = c(64, 64))
# fit <- npe(prior, simulator, density_estimator = "maf", embedding_net = emb)
```
