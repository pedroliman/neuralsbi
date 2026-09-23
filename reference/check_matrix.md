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
check_matrix(value, d = NULL, arg, what = NULL, min_rows = 0L)
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

- min_rows:

  Smallest number of rows accepted. `0` (the default) accepts an empty
  matrix, which is right for most callers –
  [`within_support()`](https://neuralsbi.pedrodelima.com/reference/within_support.md)
  on zero draws is a legitimate, if odd, question. An observation
  argument is not:
  [`resolve_obs()`](https://neuralsbi.pedrodelima.com/reference/resolve_obs.md)
  passes `min_rows = 1` so a zero-row `obs`/`x` (a real-data filter that
  happened to drop every row) gets a named error here instead of
  reaching `x[1, , drop = FALSE]` and failing with a bare, unnamed
  "subscript out of bounds" (#349).

## Value

A numeric matrix with `d` columns, column names preserved.
