# Reject sampler arguments that `sampler` does not read

[`mcmc_posterior()`](https://neuralsbi.pedrodelima.com/reference/mcmc_posterior.md)
keeps `...` as a list and the samplers read it by name, so a misspelled
argument (`wamrup`) was dropped and the default used with no message
(#390).

## Usage

``` r
check_sampler_dots(dots, sampler)
```

## Arguments

- dots:

  The sampler arguments collected from `...`.

- sampler:

  Either `"slice"` or `"stan"`.
