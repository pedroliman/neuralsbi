# Validate a callback argument

A function stored now and called later fails at the call site, which can
be a training run and several frames away from the argument that was
wrong. Checking at the constructor that it is a function, and that it
can take the argument it will be given, keeps the complaint next to the
mistake.

## Usage

``` r
check_function(f, arg, what = NULL)
```

## Arguments

- f:

  The user's value.

- arg:

  Name of the argument.

- what:

  Optional phrase naming the argument the function receives, e.g.
  `"the number of draws"`. Shown in parentheses.

## Value

`f`, invisibly.
