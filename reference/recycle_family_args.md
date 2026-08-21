# Recycle a family's arguments to one value per parameter

`prior_gamma(shape = c(2, 3), rate = 1)` is two parameters with a shared
rate, and `prior_gamma(shape = c(2, 3), rate = c(1, 2, 3))` is a
mistake. The longest argument sets the number of parameters and every
other one has to be that length or length 1, which is
[`prior_normal()`](https://neuralsbi.pedrodelima.com/reference/prior_normal.md)'s
rule for `sd` generalized to families with more than one parameter.

## Usage

``` r
recycle_family_args(args, fname)
```

## Arguments

- args:

  Named list of numeric vectors.

- fname:

  Constructor name, for the error message.

## Value

`args`, each entry a length-`d` double vector without names.
