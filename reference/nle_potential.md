# Unnormalized log posterior of an NLE fit

\\\log q\_\phi(x \mid \theta) + \log p(\theta)\\, returning `-Inf`
outside the prior support. This is the potential the MCMC samplers
target.

## Usage

``` r
nle_potential(fit, x_obs, max_batch = 1e+05)
```
