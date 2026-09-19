# Validate a flag argument

One `TRUE` or `FALSE`. Branching on a flag with
[`isTRUE()`](https://rdrr.io/r/base/Logic.html) alone accepts the value
silently but only ever takes its "false" branch for anything that is not
the literal value `TRUE` – `"yes"`, `1`, or `"TRUE"` as a string all
pass through and flip what the caller gets back with no error (#321).
`check_flag()` rejects those up front so the mistake is caught at the
call that made it rather than read off a wrong-shaped result.

## Usage

``` r
check_flag(x, arg)
```

## Arguments

- x:

  The user's value.

- arg:

  Name of the argument.

## Value

`x`, unchanged.
