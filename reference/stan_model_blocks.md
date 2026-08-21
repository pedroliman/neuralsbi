# Data, parameters and model blocks around the generated likelihood

[`prior_uniform()`](https://neuralsbi.pedrodelima.com/reference/prior_uniform.md)
and
[`prior_normal()`](https://neuralsbi.pedrodelima.com/reference/prior_normal.md)
are written out through the data block, so one compiled model serves any
bounds and any prior mean. Every other named family (see
[prior_families](https://neuralsbi.pedrodelima.com/reference/prior_families.md))
becomes a literal Stan sampling statement instead: its parameters are
what pick the distribution, and there is nothing left to vary at run
time. A prior built with
[`prior_custom()`](https://neuralsbi.pedrodelima.com/reference/prior_custom.md)
is arbitrary R code with no Stan counterpart, so the model block is the
user's to write.

## Usage

``` r
stan_model_blocks(fit, name, packed)
```
