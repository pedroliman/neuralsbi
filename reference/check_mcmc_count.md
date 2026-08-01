# Validate one MCMC count argument and return it as an integer

[`as.integer()`](https://rdrr.io/r/base/integer.html) on its own lets a
nonsensical count through. `thin = 0` is the one that bites:
[`slice_sample()`](https://neuralsbi.pedrodelima.com/reference/slice_sample.md)
then runs `warmup` iterations and no kept ones, and returns its
zero-initialized array of draws with diagnostics computed on it. The
counts are checked here, before anything is stored on the posterior
object.

## Usage

``` r
check_mcmc_count(value, name, min, why = NULL)
```

## Arguments

- value:

  The supplied value.

- name:

  Argument name, for the error message.

- min:

  Smallest value allowed.

- why:

  Optional clause explaining the bound.
