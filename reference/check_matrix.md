# Validate a matrix argument at a public boundary

Unlike
[`as_theta_matrix()`](https://neuralsbi.pedrodelima.com/reference/as_theta_matrix.md),
which reshapes whatever it is given, this errors rather than guess. A
bare vector is read as a single row and must therefore have exactly `d`
entries; a length that does not match used to be recycled into a matrix
of the right width, which turned one wrong-length parameter vector into
several parameter sets and a plausible-looking answer.

## Usage

``` r
check_matrix(value, d = NULL, arg, what = NULL)
```

## Arguments

- value:

  The user's value: a numeric vector, matrix or data frame.

- d:

  Required number of columns, or `NULL` to accept any width.

- arg:

  Name of the argument, as it appears in the user's call.

- what:

  Optional phrase describing what a column means, e.g.
  `"one parameter per column"`. Shown in parentheses.

## Value

A numeric matrix with `d` columns, column names preserved.
