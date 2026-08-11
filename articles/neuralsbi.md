# Getting started with neuralsbi

This tutorial provides a short introduction to neural simulation-based
inference (SBI) in R using the `neuralsbi` R package. Here we will focus
on one of three key sets of methods in neural SBI: neural posterior
estimation (NPE). The goal in NPE is to approximate the posterior
distribution p(\theta \| x ) of parameters of interest \theta. Here we
will use the example of stochastic SIR model.

## Inputs to neural posterior estimation: A simulator and a prior distribution

`neuralsbi` requires two inputs: a prior over the parameters you want to
estimate, and a simulator that maps one parameter vector set into one
vector of outcomes. The simulator can (and normally should) be
stochastic, so the same parameter vector \theta can result in many
possible realizations. Variation in observed outcomes can come from both
the underlying process (i.e., disease spreads stochastically) and
measurement (i.e., we don’t observe all cases).

### Define your simulator

Here we write a “chain binomial” SIR model in daily time-steps, with a
population of size N. The parameters we want to estimate are \beta, the
recovery rate \gamma, and we see a fraction \rho of the new infections
as reported cases. In reality, we would want to inform \gamma using
external information, which we can do by providing an informative prior
for \gamma.

``` r

library(neuralsbi)
#> 
#> Attaching package: 'neuralsbi'
#> The following object is masked from 'package:base':
#> 
#>     sample
library(future)

# TODO: Just write the model with a daily time-step, don't do this week and daily loops. Just aggregate weekly observations after the fect.
# Also, don't return log1p as cases reported, we shoulnd't have to do that.
# # Also, let's just return the daily reported cases.
# TODO: Let's replace the weeks argument with "days"
sir_simulator <- function(beta, gamma, rho, N = 1e5, I0 = 20, weeks = 8) {
  S <- N - I0; I <- I0
  reported <- numeric(weeks)
  for (w in seq_len(weeks)) {
    infections <- 0
    for (d in seq_len(7)) {                             # a day at a time
      new_inf <- rbinom(1, S, 1 - exp(-beta * I / N))
      new_rec <- rbinom(1, I, 1 - exp(-gamma))
      S <- S - new_inf
      I <- I + new_inf - new_rec
      infections <- infections + new_inf
    }
    reported[w] <- rbinom(1, infections, rho)           # only rho of them are seen
  }
  # TODO
  log1p(reported)
}

# Define our priors:
prior <- prior_uniform(low  = c(beta = 0.20, gamma = 0.08, rho = 0.1),
                       high = c(beta = 0.60, gamma = 0.20, rho = 0.7))
```

It is always a good idea to test our mode, so here we simulate an
outbreak to test our model:
