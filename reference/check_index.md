# Validate a column index, which may be given as a name

`plot_sbc(sbc_result, param = 99)` used to reach `ranks[, 99]` and
report "subscript out of bounds", which names neither the argument nor
how many parameters there are. The rank matrix carries `colnames` and
every other plotting function labels by name, so a name is accepted here
too and resolved against `nms`.

## Usage

``` r
check_index(value, arg, nms = NULL, n, what = "parameter")
```

## Arguments

- value:

  The user's value: one index, or one name to match against `nms`.

- arg:

  Name of the argument.

- nms:

  Names to match a character value against, or `NULL` when the columns
  are unnamed.

- n:

  Number of columns.

- what:

  Noun for one column, e.g. `"parameter"`.

## Value

An integer index between 1 and `n`.
