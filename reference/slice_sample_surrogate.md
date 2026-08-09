# Slice-sample the unnormalized posterior of an [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) or [`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md) fit

The sampler does not care which surrogate produced the potential, so
[`surrogate_potential()`](https://neuralsbi.pedrodelima.com/reference/surrogate_potential.md)
is the only line that looks at which fit it has.

## Usage

``` r
slice_sample_surrogate(fit, x_obs, ctl, n, verbose = FALSE)
```
