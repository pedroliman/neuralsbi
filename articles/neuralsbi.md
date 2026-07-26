# Getting started with neuralsbi

`neuralsbi` performs **Neural Posterior Estimation (NPE)**: given a
prior and a simulator, it trains a neural network to approximate the
Bayesian posterior `p(theta | x)` — no likelihood required. This
vignette walks through the core workflow and shows how to check that the
result is trustworthy.

## The three ingredients

Every SBI problem needs (1) a **prior**, (2) a **simulator**, and (3) an
**observation**. Here is a small linear-Gaussian model whose posterior
we happen to know in closed form, so we can check our answer.

``` r

library(neuralsbi)
#> 
#> Attaching package: 'neuralsbi'
#> The following object is masked from 'package:base':
#> 
#>     sample

# (1) prior over two parameters
prior <- prior_normal(mean = c(mu1 = 0, mu2 = 0), sd = 1)

# (2) simulator: one parameter set in, one observation out. The prior names
# match the arguments, so mu1 and mu2 arrive by name.
sigma <- 0.5
simulator <- function(mu1, mu2) c(mu1, mu2) + rnorm(2, sd = sigma)

# (3) the observation we want to explain
x_obs <- c(1.0, -0.5)
```

## Train an amortized posterior

[`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) draws
parameters from the prior, runs the simulator, and trains a conditional
density estimator. For this linear-Gaussian model, we use the
closed-form conditional-Gaussian estimator, which is *exact* and
requires no neural network training.

``` r

fit <- npe(prior, simulator, n_simulations = 2000,
           density_estimator = "linear_gaussian", seed = 1)
fit
#> <nsbi_npe> Neural Posterior Estimation fit
#>   density estimator : linear_gaussian
#>   parameters (dim)  : 2
#>     names           : mu1, mu2 
#>   data (dim)        : 2
#>   simulations       : 2000
#>   -> build a posterior with posterior(fit, x_obs = ...)
```

The fitted network approximates the posterior for *any* observation, not
just the one at hand — this is what “amortized” means. You train once,
then condition on new data without refitting.

``` r

post  <- posterior(fit, x_obs = x_obs)
draws <- sample(post, 2000)

colMeans(draws)          # posterior mean
#>        mu1        mu2 
#>  0.7951043 -0.4207836
map_estimate(post)       # MAP point estimate
#>        mu1        mu2 
#>  0.8023301 -0.4156778
pairplot(draws)          # joint + marginal view
```

![Pairs plot of posterior draws for the two
parameters.](figures/neuralsbi-unnamed-chunk-4-1.png)

plot of chunk unnamed-chunk-4

## Did it work? Check against the truth

For this model the posterior is a known Gaussian. Let us compare.

``` r

d <- 2
Sigma <- solve(diag(d) + diag(d) / sigma^2)
mu    <- as.numeric(Sigma %*% (x_obs / sigma^2))

rbind(analytic = mu, estimated = colMeans(draws))
#>                 mu1        mu2
#> analytic  0.8000000 -0.4000000
#> estimated 0.7951043 -0.4207836

# classifier two-sample test: ~0.5 => our samples look like analytic samples
z <- matrix(rnorm(2000 * d), ncol = d)
analytic_draws <- sweep(z %*% chol(Sigma), 2, mu, `+`)
c2st(draws, analytic_draws)$accuracy
#> [1] 0.5055
```

## Calibration when you *don’t* know the truth

Usually there is no analytic posterior. **Simulation-based calibration
(SBC)** still tells you whether the posterior is well calibrated: it
should produce uniform rank statistics.

``` r

res <- sbc(fit, simulator, n_sbc = 100, n_posterior_samples = 100, seed = 2)
expected_coverage(res)        # nominal vs empirical credible-interval coverage
#>    nominal  mu1  mu2
#> 1     0.05 0.05 0.03
#> 2     0.10 0.09 0.06
#> 3     0.15 0.15 0.12
#> 4     0.20 0.16 0.20
#> 5     0.25 0.21 0.23
#> 6     0.30 0.26 0.28
#> 7     0.35 0.30 0.34
#> 8     0.40 0.34 0.37
#> 9     0.45 0.46 0.44
#> 10    0.50 0.50 0.47
#> 11    0.55 0.55 0.51
#> 12    0.60 0.61 0.56
#> 13    0.65 0.66 0.60
#> 14    0.70 0.75 0.65
#> 15    0.75 0.82 0.71
#> 16    0.80 0.85 0.75
#> 17    0.85 0.88 0.77
#> 18    0.90 0.94 0.86
#> 19    0.95 0.98 0.95
plot_sbc(res)                 # flat histogram = calibrated
```

![SBC rank histogram; a flat histogram indicates
calibration.](figures/neuralsbi-unnamed-chunk-6-1.png)

plot of chunk unnamed-chunk-6

## Neural estimators with torch

For non-Gaussian posteriors, you need a neural density estimator like
the **Mixture Density Network (MDN)**. This requires `torch`.

``` r

fit_mdn <- npe(prior, simulator, n_simulations = 2000,
               density_estimator = "mdn", max_epochs = 200, seed = 1)
fit_mdn
#> <nsbi_npe> Neural Posterior Estimation fit
#>   density estimator : mdn
#>   parameters (dim)  : 2
#>     names           : mu1, mu2 
#>   data (dim)        : 2
#>   simulations       : 2000
#>   best val loss     : 1.2096
#>   -> build a posterior with posterior(fit, x_obs = ...)

# on this Gaussian model the MDN matches the exact linear_gaussian posterior
draws_mdn <- sample(posterior(fit_mdn, x_obs = x_obs), 2000)
rbind(analytic = mu, mdn = colMeans(draws_mdn))
#>                mu1        mu2
#> analytic 0.8000000 -0.4000000
#> mdn      0.8393949 -0.4265376
c2st(draws_mdn, analytic_draws)$accuracy
#> [1] 0.525
```

## Non-Gaussian posteriors

The MDN is not limited to Gaussian posteriors: it recovers the bimodal,
crescent-shaped posterior of the classic **two-moons** task, and the
flow estimators (`"maf"`, `"nsf"`) go further still. That comparison is
the subject of
[`vignette("density-estimators")`](https://neuralsbi.pedrodelima.com/articles/density-estimators.md).

## Where to go next

The vignettes build on each other:

1.  [`vignette("density-estimators")`](https://neuralsbi.pedrodelima.com/articles/density-estimators.md)
    — which estimator to use, and when.
2.  [`vignette("diagnostics")`](https://neuralsbi.pedrodelima.com/articles/diagnostics.md)
    — calibration and predictive checks for a fitted posterior.
3.  [`vignette("sir-epidemic")`](https://neuralsbi.pedrodelima.com/articles/sir-epidemic.md)
    — the complete Bayesian workflow on an applied epidemic-model
    calibration.

[`?npe`](https://neuralsbi.pedrodelima.com/reference/npe.md),
[`?posterior`](https://neuralsbi.pedrodelima.com/reference/posterior.md),
and [`?sbc`](https://neuralsbi.pedrodelima.com/reference/sbc.md)
document every argument.
