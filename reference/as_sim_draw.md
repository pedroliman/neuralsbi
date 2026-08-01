# Coerce one simulator return value to a named numeric vector

Every rejection names the actual problem, because the alternative is a
coercion error from three frames down inside
[`as.matrix()`](https://rdrr.io/r/base/matrix.html).

## Usage

``` r
as_sim_draw(out, i = 1L)
```

## Arguments

- out:

  What the simulator returned.

- i:

  Index of the simulation, for the error message.
