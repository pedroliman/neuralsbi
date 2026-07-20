% Hand-written (mirrors roxygen tags); regenerate with devtools::document()
\name{tasks}
\alias{tasks}
\alias{task_gaussian_linear}
\alias{task_two_moons}
\alias{task_slcp}
\alias{task_sir}
\title{Benchmark tasks}
\usage{
task_gaussian_linear(dim = 10L, prior_var = 0.1, noise_var = 0.1)

task_two_moons()

task_slcp()

task_sir(N = 1e6, days = 160, n_points = 10L, n_obs_draws = 1000L)
}
\arguments{
\item{dim}{Parameter/data dimension for the Gaussian linear task (sbibm
uses 10).}

\item{prior_var, noise_var}{Prior and likelihood variances for the Gaussian
linear task.}

\item{N, days, n_points, n_obs_draws}{SIR task: population size, horizon in
days, number of observation times, and binomial trials per observation.}
}
\value{
A list of class \code{nsbi_task} with elements \code{name}, \code{prior},
\code{simulator}, \code{dim_theta}, \code{dim_x}, and optionally
\code{reference_posterior(x_obs, n)} returning exact posterior draws.
}
\description{
Standard simulation-based-inference benchmark tasks, matching the
definitions in the \code{sbibm} suite that the Python \code{sbi} package
benchmarks against. Each task bundles a prior, a simulator, and (where one
exists) an analytic reference posterior, so the same object drives unit
tests, calibration studies, and the head-to-head comparison with Python
\code{sbi} in \code{inst/benchmarks/}.
\itemize{
\item \code{task_gaussian_linear()} -- conjugate Gaussian; analytic posterior.
\item \code{task_two_moons()} -- crescent-shaped, bimodal posterior.
\item \code{task_slcp()} -- 5 parameters, 8-dimensional data, strongly
  non-Gaussian posterior.
\item \code{task_sir()} -- SIR epidemic dynamics with log-normal priors;
  no closed-form posterior (verify with SBC / coverage).
}
}
