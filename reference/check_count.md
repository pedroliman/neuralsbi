# Validate a count argument

One finite whole number, at least `min`.
[`as.integer()`](https://rdrr.io/r/base/integer.html) on its own accepts
`2.7` and `TRUE` and turns `NA` into a silent shortfall, so counts that
decide how many simulations to run or how many draws to keep go through
here instead.

## Usage

``` r
check_count(n, arg, min = 1L, why = NULL)
```

## Arguments

- n:

  The user's value.

- arg:

  Name of the argument.

- min:

  Smallest allowed value.

- why:

  Optional clause explaining the bound, appended to the message.

## Value

`n` as an integer.
