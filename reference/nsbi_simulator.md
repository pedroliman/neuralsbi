# The simulator contract

A simulator in `neuralsbi` is called **once per parameter set** and
returns **one simulated observation**. Everything that calls a simulator
– [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md),
[`simulate_for_sbi()`](https://neuralsbi.pedrodelima.com/reference/simulate_for_sbi.md),
[`npe_sequential()`](https://neuralsbi.pedrodelima.com/reference/npe_sequential.md),
[`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md),
[`tarp()`](https://neuralsbi.pedrodelima.com/reference/tarp.md),
[`posterior_predictive()`](https://neuralsbi.pedrodelima.com/reference/posterior_predictive.md)
– uses this contract.

## How parameters arrive

Two signatures are accepted, and which one applies is decided once per
run from `formals(simulator)`. If every parameter name appears among the
simulator's formals, the parameters are passed **by name**, one scalar
each:

    prior <- prior_uniform(low = c(mu = -5, sigma = 0.1),
                           high = c(mu =  5, sigma = 3))

    simulator <- function(mu, sigma) {
      y <- rnorm(100, mean = mu, sd = sigma)
      c(mean = mean(y), sd = sd(y))
    }

Otherwise the whole **named parameter vector** goes to the first
argument:

    simulator <- function(theta) {
      y <- rnorm(100, mean = theta["mu"], sd = theta["sigma"])
      c(mean = mean(y), sd = sd(y))
    }

An unnamed prior always takes the vector form, and so does a prior whose
names are not syntactic R names (`"beta[1]"`), since those can never
match a formal.

## Everything else the simulator needs

Observed data, a time grid, a population size, a design matrix, a solver
tolerance: anything that is not a calibrated parameter goes in
`sim_args`, a named list forwarded to every simulator call.

    simulator <- function(alpha, beta, sigma, x_grid) {
      rnorm(length(x_grid), mean = alpha + beta * x_grid, sd = sigma)
    }

    fit <- npe(prior, simulator, n_simulations = 10000,
               sim_args = list(x_grid = seq(-1, 1, length.out = 50)))

A list is used rather than `...` because every one of
[`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md)'s own
arguments would otherwise be exposed to R's partial matching: `x`,
`theta`, `n` and `seed` are all natural names for a simulator argument
and all of them would be captured silently.

## What the simulator returns

One observation per call: a numeric vector of length `d`, a scalar, or a
one-row matrix or data frame. Names on the vector, or `colnames` on the
matrix or data frame, become the outcome names used in summaries and
plots. Every draw must agree on `d`. Lists, character or factor columns,
and anything with more than one row are rejected by name.

One sharp edge: arithmetic on a single named parameter keeps the name,
so `theta["beta"] / theta["gamma"]` returns a scalar named `beta` and
would quietly become an outcome called `beta`. Vectors are unaffected –
`theta["alpha"] + theta["beta"] * x` comes back unnamed. Use
[`unname()`](https://rdrr.io/r/base/unname.html), or name the outcome
deliberately with `c(r0 = ...)`.

## Failed simulations

A draw whose output contains `NA`, `NaN` or an infinite value is
dropped, together with its parameters, and one warning per run reports
the count and the rate. A single non-finite value would otherwise poison
the training loss for the whole fit and surface much later as a `NaN`
validation loss.

Dropping conditions on the simulator having succeeded. When failure
depends on the parameters – an ODE that diverges for large `beta`, a
model that returns `NA` outside a stability region – the surviving draws
are no longer a sample from the prior, and the fit targets the posterior
*given success* rather than the posterior. A high drop rate is a
modelling signal, and the honest fix is usually to handle the failure
inside the simulator.

## See also

[nsbi_parallel](https://neuralsbi.pedrodelima.com/reference/nsbi_parallel.md)
for how the calls are spread across workers.
