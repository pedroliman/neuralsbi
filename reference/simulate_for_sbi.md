# Run a simulator over prior draws

Draws `n` parameter vectors from the prior and calls the simulator once
per draw. Under a future plan the draws are spread across workers. See
[nsbi_simulator](https://neuralsbi.pedrodelima.com/reference/nsbi_simulator.md),
[nsbi_parallel](https://neuralsbi.pedrodelima.com/reference/nsbi_parallel.md)
and
[nsbi_progress](https://neuralsbi.pedrodelima.com/reference/nsbi_progress.md).

## Usage

``` r
simulate_for_sbi(
  simulator,
  prior,
  n,
  sim_args = list(),
  seed = NULL,
  verbose = FALSE
)
```

## Arguments

- simulator:

  A function called once per parameter set, returning one simulated
  observation: a numeric vector, a scalar, or a one-row matrix or data
  frame. See
  [nsbi_simulator](https://neuralsbi.pedrodelima.com/reference/nsbi_simulator.md).

- prior:

  An `nsbi_prior` (see
  [`prior_uniform()`](https://neuralsbi.pedrodelima.com/reference/prior_uniform.md),
  [`prior_normal()`](https://neuralsbi.pedrodelima.com/reference/prior_normal.md)).

- n:

  Number of simulations.

- sim_args:

  Named list of extra arguments passed to every simulator call: observed
  data, a time grid, a fixed population size, solver settings. See
  [nsbi_simulator](https://neuralsbi.pedrodelima.com/reference/nsbi_simulator.md).

- seed:

  Optional integer seed for reproducibility.

- verbose:

  Print training progress.

## Value

A list with `theta` (`n x dim`) and `x` (`n x d`) matrices and
`n_dropped`, the number of simulations discarded for non-finite output.

## Details

Simulations whose output is not finite are dropped together with their
parameters, with a warning.

## Examples

``` r
prior <- prior_uniform(c(a = -1, b = -1), c(a = 1, b = 1))
sims <- simulate_for_sbi(function(a, b) c(a^2, b^2), prior, n = 100)
str(sims)
#> List of 3
#>  $ theta    : num [1:100, 1:2] 0.0598 0.2033 -0.1886 0.8333 -0.9285 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:2] "a" "b"
#>  $ x        : num [1:100, 1:2] 0.00358 0.04134 0.03558 0.69433 0.86214 ...
#>  $ n_dropped: int 0
```
