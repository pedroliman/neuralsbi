# Drop simulations whose parameters or output are not finite

A row survives only when every value in both `theta` and `x` is finite,
and the two are subset together so the training set stays paired.
Checking `theta` matters on the pre-computed path, where a user-supplied
`NA` would otherwise reach [`chol()`](https://rdrr.io/r/base/chol.html)
in the estimator and surface as a linear-algebra error. Returns the
surviving pair and the number dropped; warns once, and errors when
nothing survives. See the "Failed simulations" section of
[nsbi_simulator](https://neuralsbi.pedrodelima.com/reference/nsbi_simulator.md).

## Usage

``` r
drop_failed_sims(
  theta,
  x,
  what = "simulations",
  x_hint = "Check the simulator on a single prior draw."
)
```

## Arguments

- theta:

  Parameter matrix, or `NULL` when the caller has no parameters to keep
  aligned.

- x:

  Simulated data matrix.

- what:

  Noun for the warning ("simulations", "SBC trials", ...).

- x_hint:

  What to check when the non-finite values are in `x`. The default
  points at the simulator, which is only useful advice when one was run;
  the pre-computed `(theta, x)` path passes advice about `x` itself.
