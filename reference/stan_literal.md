# One prior parameter or bound as a Stan literal

Not `stan_num()`: that one writes `%.17g` because a trained weight has
to come back bit for bit, and `0.40000000000000002` in the middle of a
`lognormal(...)` call is noise in a block a user is meant to read.
Fifteen significant digits is more precision than any hand-written prior
parameter carries.

## Usage

``` r
stan_literal(x)
```

## Arguments

- x:

  One number.
