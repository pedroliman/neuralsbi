# Save and reload a fitted model

A fit whose estimator is `"maf"`, `"mdn"` or `"nsf"` holds a `torch`
module, and a torch module is an external pointer.
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) writes the pointer,
not the network: the file reloads without complaint and the object
prints normally, but the first call that touches the network fails with
`external pointer is not valid`. `save_npe()` and `load_npe()` are the
round trip that works.

## Usage

``` r
save_npe(fit, path)

load_npe(path)

save_nle(fit, path)

load_nle(path)
```

## Arguments

- fit:

  An `nsbi_npe` object from
  [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) or
  [`npe_sequential()`](https://neuralsbi.pedrodelima.com/reference/npe_sequential.md),
  or an `nsbi_nle` object from
  [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md).

- path:

  File to write to (or read from). The convention is `.rds`.

## Value

`save_npe()` returns `path` invisibly. `load_npe()` returns the fit.

## Details

`save_npe()` writes one `.rds` file holding the network's weights (via
[`torch::torch_save()`](https://torch.mlverse.org/docs/reference/torch_save.html)
on its `state_dict`) alongside everything else the fit carries as
ordinary R objects: the prior, the standardization centers and scales,
parameter and outcome names, the simulation count, and the training
history. `load_npe()` rebuilds the network from the recorded
architecture and restores the weights, returning an `nsbi_npe` that
behaves exactly like the one you trained.

A `"linear_gaussian"` fit holds no torch objects and round-trips through
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) unharmed;
`save_npe()` accepts it anyway, so saving code does not have to know
which estimator was used.

Weights are saved, not code. A fit saved by one version of `neuralsbi`
loads into a later one as long as the estimator's architecture has not
changed; `load_npe()` reports the version that wrote the file when the
rebuild fails.

`save_nle()` and `load_nle()` are aliases; both pairs handle either kind
of fit, and the two names exist only so calling code reads the way the
fit was made.

## Examples

``` r
prior <- prior_uniform(c(mu = -2, nu = -2), c(mu = 2, nu = 2))
simulator <- function(mu, nu) c(a = mu + rnorm(1, sd = 0.1),
                                b = nu + rnorm(1, sd = 0.1))
fit <- npe(prior, simulator, n_simulations = 500,
           density_estimator = "linear_gaussian")

path <- tempfile(fileext = ".rds")
save_npe(fit, path)
fit2 <- load_npe(path)
sample(posterior(fit2, x_obs = c(0.8, 0.6)), 100)
#> <nsbi_samples> 100 draws x 2 parameters
#>   support acceptance rate: 1.000
#>             mu        nu
#> [1,] 0.7834703 0.6850203
#> [2,] 0.9073065 0.5592285
#> [3,] 0.7287850 0.5026705
#> [4,] 0.6951093 0.3997085
#> [5,] 0.7223549 0.6151154
#> [6,] 0.7333032 0.6719495
unlink(path)
```
