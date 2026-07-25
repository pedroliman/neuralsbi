# Run a simulator over prior draws

Draws `n` parameter vectors from the prior and pushes them through the
simulator. The simulator is called on chunks of rows rather than on the
whole matrix at once, which is what lets the run report progress and,
under a future plan, spread the chunks across workers. See
[nsbi_parallel](https://pedroliman.github.io/neuralsbi/reference/nsbi_parallel.md)
and
[nsbi_progress](https://pedroliman.github.io/neuralsbi/reference/nsbi_progress.md).

## Usage

``` r
simulate_for_sbi(
  simulator,
  prior,
  n,
  seed = NULL,
  verbose = FALSE,
  chunk_size = NULL
)
```

## Arguments

- simulator:

  A function mapping an `n x dim` matrix of parameters to an `n x d`
  matrix of simulated data. Ignored if `theta` and `x` are given. Column
  names on its output (e.g. via `colnames(out) <- c("cases_wk1", ...)`)
  become the outcome names used in plots.

- prior:

  An `nsbi_prior` (see
  [`prior_uniform()`](https://pedroliman.github.io/neuralsbi/reference/prior_uniform.md),
  [`prior_normal()`](https://pedroliman.github.io/neuralsbi/reference/prior_normal.md)).

- n:

  Number of simulations.

- seed:

  Optional integer seed for reproducibility.

- verbose:

  Print training progress.

- chunk_size:

  Rows of `theta` per simulator call. `NULL` (default) splits the run
  into about 64 chunks. Chunks are the unit of work sent to `future`
  workers and the unit of progress reporting; see
  [nsbi_parallel](https://pedroliman.github.io/neuralsbi/reference/nsbi_parallel.md).

## Value

A list with `theta` (`n x dim`) and `x` (`n x d`) matrices.

## Examples

``` r
prior <- prior_uniform(c(-1, -1), c(1, 1))
sims <- simulate_for_sbi(function(theta) theta^2, prior, n = 100)
str(sims)
#> List of 2
#>  $ theta: num [1:100, 1:2] -0.0435 0.4745 -0.8025 -0.3234 0.3867 ...
#>  $ x    : num [1:100, 1:2] 0.00189 0.22512 0.64396 0.10456 0.14956 ...
```
