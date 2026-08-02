# Call the simulator for one parameter set, naming it if it fails

[`as_sim_draw()`](https://neuralsbi.pedrodelima.com/reference/as_sim_draw.md)
reports everything wrong with a simulator's *return value* and names the
simulation it came from. A failure inside the simulator happens before
there is a return value, so on its own it arrives as bare as R left it:
`argument "ls" is missing, with no default`, with no index, no
parameters, and nothing to suggest that the real problem is a dispatch
mismatch rather than a missing default. Under a future plan the call
crosses a worker boundary first, which is where a useful traceback goes.

## Usage

``` r
call_sim_once(call_one, theta_i, i = 1L)
```

## Arguments

- call_one:

  The closure from
  [`simulator_caller()`](https://neuralsbi.pedrodelima.com/reference/simulator_caller.md).

- theta_i:

  One parameter set.

- i:

  Index of the simulation, for the error message.

## Value

Whatever the simulator returned.

## Details

The original condition is kept as the parent, and its class is carried
onto the re-raised error, so a caller catching a condition the simulator
signals still catches it.
