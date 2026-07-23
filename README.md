
<!-- README.md is generated from README.Rmd. Please edit that file, then run rmarkdown::render("README.Rmd"). -->

# neuralsbi

<!-- badges: start -->

[![R-CMD-check](https://github.com/pedroliman/neuralsbi/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pedroliman/neuralsbi/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/pedroliman/neuralsbi/actions/workflows/pkgdown.yaml/badge.svg)](https://pedroliman.github.io/neuralsbi/)
<!-- badges: end -->

`neuralsbi` is an R-native package for [Neural Simulation-based
inference](https://simulation-based-inference.org).

Neural estimators are implemented directly in R on the
[`torch`](https://torch.mlverse.org/) R package.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("pedroliman/neuralsbi")

# the neural back end (once)
install.packages("torch")
torch::install_torch()
```

## Usage

Simulation-based inference fits a posterior from a prior and a
simulator, with no likelihood required. To keep the setup familiar, here
is ordinary linear regression written as a simulator: a response `y`
scattered around a line, `y ~ Normal(alpha + beta * x, sigma)` — the
same model you might write in Stan. Its likelihood is easy to write
down, which is exactly what makes it a good check: we know the
coefficients that generated the data, so we can confirm the posterior
recovers them.

``` r
library(neuralsbi)
set.seed(1)

# Regression design: a single covariate x, measured at 50 fixed points.
N <- 50
x <- seq(-1, 1, length.out = N)

# Simulator: given rows of (alpha, beta, sigma), draw one response vector y each
# from y ~ Normal(alpha + beta * x, sigma). Fully vectorised over the rows of
# theta, and it only generates data — no fitting happens here.
simulator <- function(theta) {
  alpha <- theta[, 1]
  beta  <- theta[, 2]
  sigma <- theta[, 3]
  mu <- outer(beta, x) + alpha                       # row i is the line alpha_i + beta_i * x
  mu + matrix(rnorm(length(mu)), nrow(mu)) * sigma   # add row-specific Gaussian noise
}

# Priors over the intercept, slope, and noise scale, then train the posterior.
prior <- prior_uniform(low = c(-3, -3, 0.1), high = c(3, 3, 2))
fit   <- npe(prior, simulator, n_simulations = 10000, seed = 1)

# Simulate one data set from known coefficients, then infer them back. The
# observation the posterior conditions on is the response vector y.
theta_true <- c(alpha = 2, beta = -1, sigma = 0.5)
set.seed(38)                       # a fixed, representative data set
y_obs      <- simulator(rbind(theta_true))
post       <- posterior(fit, x_obs = y_obs)
draws      <- sample(post, 10000)
```

The posterior mean recovers the coefficients that generated the data:

``` r
rbind(truth = theta_true, posterior_mean = colMeans(draws))
#>                   alpha      beta     sigma
#> truth          2.000000 -1.000000 0.5000000
#> posterior_mean 2.003668 -1.041009 0.5003288
```

``` r
pairplot(draws, truth = theta_true, labels = c("alpha", "beta", "sigma"))
```

<img src="man/figures/README-pairplot-1.png" alt="Pairwise posterior over the regression intercept, slope, and noise scale." width="100%" />

The same posterior gives a point estimate; calibration checks such as
simulation-based calibration live in `vignette("diagnostics")`.

``` r
map_estimate(post)     # posterior mode
#> [1]  1.9989619 -1.0442036  0.4702335
```

If you’re interested in sbi in other languages or functionality not
available here, see the [awesome neural SBI
repo](https://github.com/smsharma/awesome-neural-sbi); there are some
good implementations in python and in Julia.

## Learn more

The [package website](https://pedroliman.github.io/neuralsbi/) has four
vignettes that build on each other:

1.  [Getting
    started](https://pedroliman.github.io/neuralsbi/articles/neuralsbi.html)
    — the core prior/simulator/posterior workflow.
2.  [Choosing a density
    estimator](https://pedroliman.github.io/neuralsbi/articles/density-estimators.html)
    — MDN, MAF, NSF, and the torch-free baseline.
3.  [Checking the
    posterior](https://pedroliman.github.io/neuralsbi/articles/diagnostics.html)
    — calibration and predictive diagnostics.
4.  [Case study: an SIR epidemic
    model](https://pedroliman.github.io/neuralsbi/articles/sir-epidemic.html)
    — the full Bayesian workflow on an applied problem.

## License

MIT
