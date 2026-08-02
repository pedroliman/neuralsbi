# Validate a strictly positive scalar

Validate a strictly positive scalar

## Usage

``` r
check_positive(x, arg, allow_inf = FALSE)
```

## Arguments

- x:

  The user's value.

- arg:

  Name of the argument.

- allow_inf:

  Accept `Inf`, for a bound that is disabled by setting it to infinity
  (`clip_grad_norm`).

## Value

`x` as a double.
