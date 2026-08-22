# Check that torch is available, error otherwise

Check that torch is available, error otherwise

## Usage

``` r
require_torch(
  what = "This density estimator",
  alternative = paste("Alternatively use", "density_estimator =",
    "'linear_gaussian' for a", "torch-free baseline.")
)
```

## Arguments

- what:

  What needs torch, as the subject of the sentence. The default suits
  the density estimators, which is where most of the calls are.

- alternative:

  The torch-free thing to do instead, as a full sentence. Every caller
  has one, and naming it is the difference between a dead end and a next
  step.
