# Everything both [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) and [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) do before touching a density estimator

Get the simulations (running the simulator or taking pre-computed ones),
coerce them to matrices, drop non-finite draws, and learn the two
standardizers. Shared so the two entry points cannot drift apart – they
differ only in which side of `(theta, x)` the estimator treats as its
target.

## Usage

``` r
prepare_simulations(
  prior,
  simulator,
  n_simulations,
  sim_args,
  theta,
  x,
  standardize,
  seed,
  verbose
)
```

## Value

A list with the cleaned `theta`/`x`, their standardized versions, the
two `nsbi_standardizer`s, the column names, and `n_dropped`.
