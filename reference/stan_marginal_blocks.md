# Stan blocks for a prior made of named marginals

Bounds shared by every parameter go on the `vector` declaration, which
is how the program would be written by hand. Bounds that differ cannot:
`vector<lower=...>` takes one bound for the whole container, so those
parameters are declared one at a time and assembled into `theta`
afterwards. Either way the sampling statements name a scalar, which is
what Stan's `T[,]` truncation syntax requires.

## Usage

``` r
stan_marginal_blocks(marginals, Q)
```

## Arguments

- marginals:

  The prior's marginals.

- Q:

  Number of parameters.
