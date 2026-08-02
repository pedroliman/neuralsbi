# Require numeric data, naming any column that is not

The type half of
[`check_matrix()`](https://neuralsbi.pedrodelima.com/reference/check_matrix.md),
split out because the shape rules differ between entry points but this
rule does not. On the pre-computed `(theta, x)` path a row is one
simulation, so a bare vector is a column of values rather than
[`check_matrix()`](https://neuralsbi.pedrodelima.com/reference/check_matrix.md)'s
single row; a character or factor column is the same mistake in either
place. Left to `storage.mode(x) <- "double"` such a column becomes all
`NA`, every row is then dropped as non-finite, and the error blames a
simulator that was never called.

## Usage

``` r
check_numeric(value, arg)
```

## Arguments

- value:

  The user's value: a numeric vector, matrix or data frame.

- arg:

  Name of the argument, as it appears in the user's call.

## Value

`value`, with a data frame converted to a matrix.
