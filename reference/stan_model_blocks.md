# Data, parameters and model blocks around the generated likelihood

Only
[`prior_uniform()`](https://neuralsbi.pedrodelima.com/reference/prior_uniform.md)
and
[`prior_normal()`](https://neuralsbi.pedrodelima.com/reference/prior_normal.md)
can be written out; anything built with
[`prior_custom()`](https://neuralsbi.pedrodelima.com/reference/prior_custom.md)
is arbitrary R code with no Stan counterpart, so the model block is the
user's to write.

## Usage

``` r
stan_model_blocks(fit, name, packed)
```
