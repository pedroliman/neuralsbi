# Build the closure that calls the simulator for one parameter set

The dispatch decision and the `sim_args` checks happen here, once, so
the per-draw loop is a single
[`do.call()`](https://rdrr.io/r/base/do.call.html).

## Usage

``` r
simulator_caller(simulator, param_names, sim_args = list())
```
