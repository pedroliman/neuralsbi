# Validate a vector of counts

[`check_count()`](https://neuralsbi.pedrodelima.com/reference/check_count.md)
for an argument that is a vector of counts rather than one, such as
`hidden`, where each entry is a layer width and a single bad entry
breaks the network the same way a bad scalar would. The message lists
the values rather than describing the vector, since which entry is wrong
is the thing worth reading.

## Usage

``` r
check_counts(n, arg, min = 1L, what = NULL)
```

## Arguments

- n:

  The user's value.

- arg:

  Name of the argument.

- min:

  Smallest allowed value.

- what:

  Optional phrase describing what one entry means, e.g.
  `"one hidden-layer width per entry"`. Shown in parentheses.

## Value

`n` as an integer vector.
