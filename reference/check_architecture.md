# Validate the architecture arguments shared by [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) and [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md)

All of them are checked whichever estimator was asked for, so
`n_bins = 0` is an error even under `"linear_gaussian"`, which ignores
it. A value that cannot build a network is a mistake in the call whether
or not this run would have read it.

## Usage

``` r
check_architecture(n_components, n_transforms, hidden, n_bins, tail_bound)
```

## Arguments

- n_components, hidden:

  MDN settings: number of mixture components (default 10, as in `sbi`)
  and a vector of hidden-layer widths.

- n_transforms:

  MAF/NSF setting: number of stacked autoregressive transforms (default
  5, as in `sbi`).

- n_bins, tail_bound:

  NSF settings: number of spline bins per transform and the half-width
  of the interval the spline acts on (outside it the transform is the
  identity).
