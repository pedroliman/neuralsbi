# Build a prior from one named family

The body every constructor in this file shares: validate, recycle, turn
each parameter into a marginal, hand the lot to
[`marginal_prior()`](https://neuralsbi.pedrodelima.com/reference/marginal_prior.md).

## Usage

``` r
family_prior(family, args, fname, lower = -Inf, upper = Inf, type = family)
```

## Arguments

- family:

  Family name, a key of
  [`prior_family()`](https://neuralsbi.pedrodelima.com/reference/prior_family.md).

- args:

  Named list of the validated arguments, in any order.

- fname:

  Constructor name, for error messages.

- lower, upper:

  Truncation applied to every parameter, for the half families.

- type:

  The `type` recorded on the prior; defaults to the family name.
