# Run a simulator over a parameter matrix

The single point through which every simulator call in the package
passes: the per-parameter-set contract, the future backend, RNG streams,
and the progress bar all live here. See
[nsbi_simulator](https://neuralsbi.pedrodelima.com/reference/nsbi_simulator.md)
and
[nsbi_parallel](https://neuralsbi.pedrodelima.com/reference/nsbi_parallel.md).

## Usage

``` r
run_simulator(
  simulator,
  theta,
  sim_args = list(),
  label = "Simulating",
  d = NULL
)
```

## Arguments

- simulator:

  The user's simulator, called once per row of `theta`.

- theta:

  Parameter matrix; its `colnames` are the parameter names used to
  decide how the simulator receives its arguments.

- sim_args:

  Named list of extra arguments forwarded to every call.

- label:

  Phase name for the progress bar.

- d:

  Expected number of output columns, if known.

## Value

An `n x d` matrix, one row per row of `theta`.
