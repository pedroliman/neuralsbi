# Build `sbibm`'s C2ST classifier as a torch module

A plain `linear/relu` stack ending in one logit, the same trunk
[`mlp_layers()`](https://neuralsbi.pedrodelima.com/reference/mlp_layers.md)
builds for every other estimator in the package. The one thing worth
doing by hand is the initialization: `scikit-learn` draws both weights
and biases from Glorot uniform, `torch` uses Kaiming uniform on the
weight and a fan-in bound on the bias, and the fit has to start where
`MLPClassifier`'s would for the numbers to line up.

## Usage

``` r
c2st_mlp_module(d, hidden)
```

## Arguments

- d:

  Number of columns in the draws.

- hidden:

  Widths of the hidden layers, one entry per layer.
