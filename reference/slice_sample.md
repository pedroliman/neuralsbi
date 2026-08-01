# Vectorized univariate slice sampler

Vectorized univariate slice sampler

## Usage

``` r
slice_sample(
  log_prob_fn,
  init,
  n_draws,
  warmup = 200L,
  thin = 10L,
  width = 1,
  max_steps = 100L,
  verbose = FALSE
)
```

## Arguments

- log_prob_fn:

  Vectorized log-density: takes an `m x dim` matrix, returns `m`
  log-densities (`-Inf` allowed).

- init:

  `n_chains x dim` matrix of starting points, all with finite
  log-density.

- n_draws:

  Number of retained draws in total, across all chains.

- warmup:

  Steps discarded at the start of each chain.

- thin:

  Keep one draw in `thin`.

- width:

  Initial slice width per dimension (recycled). Adapted during warmup
  and then held fixed.

- max_steps:

  Cap on stepping-out expansions per coordinate.

- verbose:

  Report progress.

## Value

A list with `draws` (`n_draws x dim`), `chains`
(`n_kept x n_chains x dim` array) and `n_evals`.
