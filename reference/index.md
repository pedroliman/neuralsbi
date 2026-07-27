# Package index

## Priors

Define the prior distribution over simulator parameters.

- [`priors`](https://neuralsbi.pedrodelima.com/reference/priors.md) :
  Priors for neural simulation-based inference
- [`prior_uniform()`](https://neuralsbi.pedrodelima.com/reference/prior_uniform.md)
  : Box-uniform (independent uniform) prior
- [`prior_normal()`](https://neuralsbi.pedrodelima.com/reference/prior_normal.md)
  : Independent normal prior
- [`prior_custom()`](https://neuralsbi.pedrodelima.com/reference/prior_custom.md)
  : Build a prior from arbitrary sampling / density functions
- [`sample_prior()`](https://neuralsbi.pedrodelima.com/reference/sample_prior.md)
  : Draw samples from a prior
- [`within_support()`](https://neuralsbi.pedrodelima.com/reference/within_support.md)
  : Test whether parameters lie within the prior support

## Simulation and training

Generate training data and fit a neural posterior estimator.

- [`simulate_for_sbi()`](https://neuralsbi.pedrodelima.com/reference/simulate_for_sbi.md)
  : Run a simulator over prior draws
- [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) : Neural
  Posterior Estimation (NPE)
- [`npe_sequential()`](https://neuralsbi.pedrodelima.com/reference/npe_sequential.md)
  : Sequential NPE with truncated-prior proposals (TSNPE)
- [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) : Neural
  Likelihood Estimation (NLE)
- [`save_npe()`](https://neuralsbi.pedrodelima.com/reference/save_npe.md)
  [`load_npe()`](https://neuralsbi.pedrodelima.com/reference/save_npe.md)
  [`save_nle()`](https://neuralsbi.pedrodelima.com/reference/save_npe.md)
  [`load_nle()`](https://neuralsbi.pedrodelima.com/reference/save_npe.md)
  : Save and reload a fitted model
- [`density_estimator`](https://neuralsbi.pedrodelima.com/reference/density_estimator.md)
  : Conditional density estimators
- [`embedding_mlp()`](https://neuralsbi.pedrodelima.com/reference/embedding_mlp.md)
  : Embedding (summary) networks for structured observations

## Neural likelihood estimation

Learn a surrogate likelihood, sample it with MCMC, and export it to
Stan.

- [`log_lik()`](https://neuralsbi.pedrodelima.com/reference/log_lik.md)
  : Evaluate a surrogate likelihood
- [`likelihood_fn()`](https://neuralsbi.pedrodelima.com/reference/likelihood_fn.md)
  : A surrogate likelihood as a plain R function
- [`posterior(`*`<nsbi_nle>`*`)`](https://neuralsbi.pedrodelima.com/reference/posterior.nsbi_nle.md)
  : Posterior from a neural likelihood
- [`sample(`*`<nsbi_nle_posterior>`*`)`](https://neuralsbi.pedrodelima.com/reference/sample.nsbi_nle_posterior.md)
  : Sample an NLE posterior with MCMC
- [`nsbi_mcmc`](https://neuralsbi.pedrodelima.com/reference/nsbi_mcmc.md)
  : MCMC over a learned likelihood
- [`stan_code()`](https://neuralsbi.pedrodelima.com/reference/stan_export.md)
  [`write_stan_model()`](https://neuralsbi.pedrodelima.com/reference/stan_export.md)
  [`stan_data()`](https://neuralsbi.pedrodelima.com/reference/stan_export.md)
  : Export a learned likelihood to Stan

## Running the simulator

The simulator contract, parallel execution, and progress reporting.

- [`nsbi_simulator`](https://neuralsbi.pedrodelima.com/reference/nsbi_simulator.md)
  : The simulator contract
- [`nsbi_parallel`](https://neuralsbi.pedrodelima.com/reference/nsbi_parallel.md)
  : Running the simulator in parallel
- [`nsbi_progress`](https://neuralsbi.pedrodelima.com/reference/nsbi_progress.md)
  [`progress-bars`](https://neuralsbi.pedrodelima.com/reference/nsbi_progress.md)
  : Progress reporting

## Working with the posterior

Condition on data, then sample, evaluate, and summarize.

- [`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md)
  : Posterior objects
- [`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md) :
  Draw samples (S3 generic)
- [`sample(`*`<nsbi_posterior>`*`)`](https://neuralsbi.pedrodelima.com/reference/sample.nsbi_posterior.md)
  : Sample from a posterior
- [`sample_posterior()`](https://neuralsbi.pedrodelima.com/reference/sample_posterior.md)
  : Sample from a posterior (non-generic alias)
- [`log_prob()`](https://neuralsbi.pedrodelima.com/reference/log_prob.md)
  : Posterior log-density
- [`map_estimate()`](https://neuralsbi.pedrodelima.com/reference/map_estimate.md)
  : Maximum a posteriori (MAP) estimate
- [`as.data.frame(`*`<nsbi_samples>`*`)`](https://neuralsbi.pedrodelima.com/reference/summaries.md)
  [`summary(`*`<nsbi_samples>`*`)`](https://neuralsbi.pedrodelima.com/reference/summaries.md)
  [`summary(`*`<nsbi_posterior>`*`)`](https://neuralsbi.pedrodelima.com/reference/summaries.md)
  [`summary(`*`<nsbi_npe>`*`)`](https://neuralsbi.pedrodelima.com/reference/summaries.md)
  : Summaries and tidy accessors

## Diagnostics

Calibration and predictive checks for a fitted posterior.

- [`diagnostics`](https://neuralsbi.pedrodelima.com/reference/diagnostics.md)
  : Posterior diagnostics
- [`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md) :
  Simulation-Based Calibration (SBC)
- [`expected_coverage()`](https://neuralsbi.pedrodelima.com/reference/expected_coverage.md)
  : Expected coverage of central credible intervals
- [`tarp()`](https://neuralsbi.pedrodelima.com/reference/tarp.md) : TARP
  expected coverage
- [`c2st()`](https://neuralsbi.pedrodelima.com/reference/c2st.md) :
  Classifier two-sample test (C2ST)
- [`posterior_predictive()`](https://neuralsbi.pedrodelima.com/reference/posterior_predictive.md)
  : Posterior predictive draws

## Plotting

- [`pairplot()`](https://neuralsbi.pedrodelima.com/reference/pairplot.md)
  : Visualize posterior samples
- [`plot_sbc()`](https://neuralsbi.pedrodelima.com/reference/plot_sbc.md)
  : Plot an SBC rank histogram
- [`plot_coverage()`](https://neuralsbi.pedrodelima.com/reference/plot_coverage.md)
  : Plot nominal vs. empirical credible-interval coverage
- [`plot_tarp()`](https://neuralsbi.pedrodelima.com/reference/plot_tarp.md)
  : Plot TARP expected coverage
- [`plot_posterior_predictive()`](https://neuralsbi.pedrodelima.com/reference/plot_posterior_predictive.md)
  : Plot posterior predictive checks

## Benchmark tasks

Reference inference problems shared by the tests and benchmarks.

- [`task_gaussian_linear()`](https://neuralsbi.pedrodelima.com/reference/tasks.md)
  [`task_two_moons()`](https://neuralsbi.pedrodelima.com/reference/tasks.md)
  [`task_slcp()`](https://neuralsbi.pedrodelima.com/reference/tasks.md)
  [`task_sir()`](https://neuralsbi.pedrodelima.com/reference/tasks.md) :
  Benchmark tasks
