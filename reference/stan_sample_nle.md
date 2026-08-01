# Sample an NLE posterior with Stan

Writes the model, compiles it, and runs NUTS. Prefers cmdstanr and falls
back to rstan.

## Usage

``` r
stan_sample_nle(fit, x_obs, ctl, n, verbose = FALSE)
```
