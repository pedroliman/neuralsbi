# Validate a file path argument

[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and
[`writeLines()`](https://rdrr.io/r/base/writeLines.html) take a
connection as well as a path, so a value that is neither gets as far as
the file system before anything complains. On the reading side
[`readRDS()`](https://rdrr.io/r/base/readRDS.html) on a path that does
not exist raises "cannot open the connection", and the file it could not
open is named only in the accompanying warning.

## Usage

``` r
check_path(path, arg = "path", must_exist = FALSE)
```

## Arguments

- path:

  The user's value.

- arg:

  Name of the argument, as it appears in the user's call.

- must_exist:

  Require the file to exist, for a path that is read.

## Value

`path`, unchanged.
