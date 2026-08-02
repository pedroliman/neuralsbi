# Validate a prior argument, and optionally its dimension

Validate a prior argument, and optionally its dimension

## Usage

``` r
check_prior(prior, arg = "prior", dim = NULL)
```

## Arguments

- prior:

  The user's value.

- arg:

  Name of the argument.

- dim:

  Number of parameters the caller expects, or `NULL` to skip that check.

## Value

The prior, invisibly.
