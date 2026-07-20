# Benchmark tasks

Standard simulation-based-inference benchmark tasks, following the
definitions in the `sbibm` benchmark suite. Each task bundles a prior, a
simulator, and (where one exists) an analytic reference posterior, so
the same object drives unit tests, calibration studies, and the
benchmark harness in `inst/benchmarks/`.

- `task_gaussian_linear()` – conjugate Gaussian; analytic posterior.

- `task_two_moons()` – crescent-shaped, bimodal posterior.

- `task_slcp()` – 5 parameters, 8-dimensional data, strongly
  non-Gaussian posterior.

- `task_sir()` – SIR epidemic dynamics with log-normal priors; no
  closed-form posterior (verify with SBC / coverage).

## Usage

``` r
task_gaussian_linear(dim = 10L, prior_var = 0.1, noise_var = 0.1)

task_two_moons()

task_slcp()

task_sir(N = 1e6, days = 160, n_points = 10L, n_obs_draws = 1000L)
```

## Arguments

- dim:

  Parameter/data dimension for the Gaussian linear task (default 10).

- prior_var, noise_var:

  Prior and likelihood variances for the Gaussian linear task.

- N, days, n_points, n_obs_draws:

  SIR task: population size, horizon in days, number of observation
  times, and binomial trials per observation.

## Value

A list of class `nsbi_task` with elements `name`, `prior`, `simulator`,
`dim_theta`, `dim_x`, and optionally `reference_posterior(x_obs, n)`
returning exact posterior draws.
