# neuralsbi

<!-- badges: start -->
[![R-CMD-check](https://github.com/pedroliman/neural.sbi/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pedroliman/neural.sbi/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/pedroliman/neural.sbi/actions/workflows/pkgdown.yaml/badge.svg)](https://pedroliman.github.io/neural.sbi/)
<!-- badges: end -->

`neuralsbi` does Bayesian inference for models with a simulator but no
tractable likelihood. Give it a prior and a function that turns parameters
into simulated data; it trains a neural conditional density estimator
(Neural Posterior Estimation) and returns a posterior you can sample from,
evaluate, and check. Neural estimators run on the
[`torch`](https://torch.mlverse.org/) R package — no Python required — and a
closed-form linear-Gaussian estimator works without `torch` at all.

## Installation

```r
# install.packages("remotes")
remotes::install_github("pedroliman/neural.sbi")

# the neural back end (once)
install.packages("torch")
torch::install_torch()
```

## Usage

```r
library(neuralsbi)

prior <- prior_uniform(low = c(-2, -2, -2), high = c(2, 2, 2))
simulator <- function(theta) {
  theta + 1 + matrix(rnorm(length(theta), sd = 0.1), nrow = nrow(theta))
}

fit   <- npe(prior, simulator, n_simulations = 2000)
post  <- posterior(fit, x_obs = c(0.8, 0.6, 0.4))
draws <- sample(post, 10000)

pairplot(draws)                # joint and marginal views
map_estimate(post)             # point estimate
sbc(fit, simulator)            # calibration check
```

The neural estimators handle genuinely non-Gaussian posteriors — here the
bimodal two-moons task:

![Two-moons posterior](man/figures/two_moons_posterior.png)

## Learn more

The [package website](https://pedroliman.github.io/neural.sbi/) has four
vignettes that build on each other:

1. [Getting started](https://pedroliman.github.io/neural.sbi/articles/neuralsbi.html) — the core prior/simulator/posterior workflow.
2. [Choosing a density estimator](https://pedroliman.github.io/neural.sbi/articles/density-estimators.html) — MDN, MAF, NSF, and the torch-free baseline.
3. [Checking the posterior](https://pedroliman.github.io/neural.sbi/articles/diagnostics.html) — calibration and predictive diagnostics.
4. [Case study: an SIR epidemic model](https://pedroliman.github.io/neural.sbi/articles/sir-epidemic.html) — the full Bayesian workflow on an applied problem.

## License

MIT
