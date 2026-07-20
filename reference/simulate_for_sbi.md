# Run a simulator over prior draws

Run a simulator over prior draws

## Usage

``` r
simulate_for_sbi(simulator, prior, n, seed = NULL, verbose = FALSE)
```

## Arguments

- simulator:

  A function mapping an `n x dim` matrix of parameters to an `n x d`
  matrix of simulated data. Ignored if `theta` and `x` are given.

- prior:

  An `nsbi_prior` (see
  [`prior_uniform()`](https://pedroliman.github.io/neural.sbi/reference/prior_uniform.md),
  [`prior_normal()`](https://pedroliman.github.io/neural.sbi/reference/prior_normal.md)).

- n:

  Number of simulations.

- seed:

  Optional integer seed for reproducibility.

- verbose:

  Print training progress.

## Value

A list with `theta` (`n x dim`) and `x` (`n x d`) matrices.
