# neuralsbi

<!-- badges: start -->
[![R-CMD-check](https://github.com/pedroliman/neuralsbi/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pedroliman/neuralsbi/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/pedroliman/neuralsbi/actions/workflows/pkgdown.yaml/badge.svg)](https://pedroliman.github.io/neuralsbi/)
<!-- badges: end -->

`neuralsbi` is an R-native package for [Neural Simulation-based inference](https://simulation-based-inference.org).

Neural estimators are implemented directly in R on the torch [`torch`](https://torch.mlverse.org/) R package.

## Installation

```r
# install.packages("remotes")
remotes::install_github("pedroliman/neuralsbi")

# the neural back end (once)
install.packages("torch")
torch::install_torch()
```

After `torch::install_torch()`, confirm the back end actually loads before
running a fit:

```r
torch::torch_tensor(1)  # should print a tensor, not an error
```

If that errors, see the troubleshooting notes below. The `linear_gaussian`
density estimator runs without torch, so you can still use the package while
you sort out the back end.

### Troubleshooting torch on macOS

On some Macs `torch::install_torch()` succeeds but the first tensor operation
fails with a message like:

```
Error: Lantern is not loaded.
... liblantern.dylib ... Symbol not found: __ZNSt13exception_ptr...
... libtorch_cpu.dylib (built for macOS 15.0 which is newer than running OS)
... Expected in: /usr/lib/libc++.1.dylib
```

The bundled libtorch was built against a newer macOS than the one you are
running, so it references a C++ standard-library symbol your system's
`libc++.1.dylib` does not have. `install_torch(reinstall = TRUE)` downloads the
same incompatible binary, so it does not help. Two things that do:

1. **Update macOS** to the version libtorch was built for (the "built for
   macOS X" line tells you which), then restart R.
2. **Install an earlier `torch`** whose libtorch targets your macOS:

   ```r
   # install.packages("remotes")
   remotes::install_version("torch", "0.13.0")
   torch::install_torch(reinstall = TRUE)
   # restart R, then: torch::torch_tensor(1)
   ```

   Try progressively older versions until `torch_tensor(1)` prints a tensor.

`neuralsbi` now detects this failure up front and prints these steps, rather
than crashing mid-fit.

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

If you're interested in sbi in other languages or functionality not available here, see the [awesome neural SBI repo](https://github.com/smsharma/awesome-neural-sbi); there are some good implementations in python and in Julia.

## Learn more

The [package website](https://pedroliman.github.io/neuralsbi/) has four
vignettes that build on each other:

1. [Getting started](https://pedroliman.github.io/neuralsbi/articles/neuralsbi.html) — the core prior/simulator/posterior workflow.
2. [Choosing a density estimator](https://pedroliman.github.io/neuralsbi/articles/density-estimators.html) — MDN, MAF, NSF, and the torch-free baseline.
3. [Checking the posterior](https://pedroliman.github.io/neuralsbi/articles/diagnostics.html) — calibration and predictive diagnostics.
4. [Case study: an SIR epidemic model](https://pedroliman.github.io/neuralsbi/articles/sir-epidemic.html) — the full Bayesian workflow on an applied problem.

## License

MIT
