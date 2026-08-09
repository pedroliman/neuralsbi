# The summary block every MCMC posterior's `print()` method shares

Everything except the closing "here is what you can do with it" lines,
which differ: an
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) fit can be
exported to Stan and an
[`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md) fit
cannot.

## Usage

``` r
cat_mcmc_posterior(x, class)
```

## Arguments

- x:

  The posterior object.

- class:

  Class name to print in the header.
