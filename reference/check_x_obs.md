# Check the observation a sequential fit targets

TSNPE spends its whole budget near one observation, so a mis-sized
`x_obs` is not a detail. Without this check
[`as_theta_matrix()`](https://neuralsbi.pedrodelima.com/reference/as_theta_matrix.md)
folds a length-3 vector into a 3 x 1 matrix and `resolve_x()` keeps row
1, so the run truncates its proposals around a value the user never
asked for and nothing says so. Called twice: once at the top of
[`npe_sequential()`](https://neuralsbi.pedrodelima.com/reference/npe_sequential.md)
for shape, and once after round 1 with `dim_x`, which is the first
moment the simulator's output width is known.

## Usage

``` r
check_x_obs(x_obs, dim_x = NULL)
```

## Arguments

- x_obs:

  The observation passed to
  [`npe_sequential()`](https://neuralsbi.pedrodelima.com/reference/npe_sequential.md).

- dim_x:

  Simulator output width, or `NULL` to check shape only.

## Value

With `dim_x`, `x_obs` as a one-row matrix; otherwise `x_obs` unchanged,
invisibly.
