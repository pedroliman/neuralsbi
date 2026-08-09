# Build the classifier torch module: `(theta, x)` -\> one logit

Three architectures behind one module, matching `sbi`'s `"resnet"`,
`"mlp"` and `"linear"` classifiers. `sbi` normalizes between the hidden
layers of its MLP (`nn.LayerNorm` by default); this one does not,
keeping the trunk the same plain `linear/relu` stack every other
estimator in the package uses
([`mlp_layers()`](https://neuralsbi.pedrodelima.com/reference/mlp_layers.md)).

## Usage

``` r
nre_module(dim_x, dim_theta, classifier, hidden, n_blocks, embedding = NULL)
```
