# Package index

## Priors

Define the prior distribution over simulator parameters.

- [`priors`](https://pedroliman.github.io/neural.sbi/reference/priors.md)
  : Priors for neural simulation-based inference
- [`prior_uniform()`](https://pedroliman.github.io/neural.sbi/reference/prior_uniform.md)
  : Box-uniform (independent uniform) prior
- [`prior_normal()`](https://pedroliman.github.io/neural.sbi/reference/prior_normal.md)
  : Independent normal prior
- [`prior_custom()`](https://pedroliman.github.io/neural.sbi/reference/prior_custom.md)
  : Build a prior from arbitrary sampling / density functions
- [`sample_prior()`](https://pedroliman.github.io/neural.sbi/reference/sample_prior.md)
  : Draw samples from a prior
- [`within_support()`](https://pedroliman.github.io/neural.sbi/reference/within_support.md)
  : Test whether parameters lie within the prior support

## Simulation and training

Generate training data and fit a neural posterior estimator.

- [`simulate_for_sbi()`](https://pedroliman.github.io/neural.sbi/reference/simulate_for_sbi.md)
  : Run a simulator over prior draws
- [`npe()`](https://pedroliman.github.io/neural.sbi/reference/npe.md) :
  Neural Posterior Estimation (NPE)
- [`npe_sequential()`](https://pedroliman.github.io/neural.sbi/reference/npe_sequential.md)
  : Sequential NPE with truncated-prior proposals (TSNPE)
- [`density_estimator`](https://pedroliman.github.io/neural.sbi/reference/density_estimator.md)
  : Conditional density estimators
- [`embedding_mlp()`](https://pedroliman.github.io/neural.sbi/reference/embedding_mlp.md)
  : Embedding (summary) networks for structured observations

## Working with the posterior

Condition on data, then sample, evaluate, and summarize.

- [`posterior()`](https://pedroliman.github.io/neural.sbi/reference/posterior.md)
  : Posterior objects
- [`sample()`](https://pedroliman.github.io/neural.sbi/reference/sample.md)
  : Draw samples (S3 generic)
- [`sample(`*`<nsbi_posterior>`*`)`](https://pedroliman.github.io/neural.sbi/reference/sample.nsbi_posterior.md)
  : Sample from a posterior
- [`sample_posterior()`](https://pedroliman.github.io/neural.sbi/reference/sample_posterior.md)
  : Sample from a posterior (non-generic alias)
- [`log_prob()`](https://pedroliman.github.io/neural.sbi/reference/log_prob.md)
  : Posterior log-density
- [`map_estimate()`](https://pedroliman.github.io/neural.sbi/reference/map_estimate.md)
  : Maximum a posteriori (MAP) estimate
- [`summary(`*`<nsbi_samples>`*`)`](https://pedroliman.github.io/neural.sbi/reference/summaries.md)
  [`summary(`*`<nsbi_posterior>`*`)`](https://pedroliman.github.io/neural.sbi/reference/summaries.md)
  [`summary(`*`<nsbi_npe>`*`)`](https://pedroliman.github.io/neural.sbi/reference/summaries.md)
  [`as.data.frame(`*`<nsbi_samples>`*`)`](https://pedroliman.github.io/neural.sbi/reference/summaries.md)
  : Summaries and tidy accessors

## Diagnostics

Calibration and predictive checks for a fitted posterior.

- [`diagnostics`](https://pedroliman.github.io/neural.sbi/reference/diagnostics.md)
  : Posterior diagnostics
- [`sbc()`](https://pedroliman.github.io/neural.sbi/reference/sbc.md) :
  Simulation-Based Calibration (SBC)
- [`expected_coverage()`](https://pedroliman.github.io/neural.sbi/reference/expected_coverage.md)
  : Expected coverage of central credible intervals
- [`tarp()`](https://pedroliman.github.io/neural.sbi/reference/tarp.md)
  : TARP expected coverage
- [`c2st()`](https://pedroliman.github.io/neural.sbi/reference/c2st.md)
  : Classifier two-sample test (C2ST)
- [`posterior_predictive()`](https://pedroliman.github.io/neural.sbi/reference/posterior_predictive.md)
  : Posterior predictive draws

## Plotting

- [`pairplot()`](https://pedroliman.github.io/neural.sbi/reference/pairplot.md)
  : Visualize posterior samples
- [`plot_sbc()`](https://pedroliman.github.io/neural.sbi/reference/plot_sbc.md)
  : Plot an SBC rank histogram
- [`plot_coverage()`](https://pedroliman.github.io/neural.sbi/reference/plot_coverage.md)
  : Plot nominal vs. empirical credible-interval coverage
- [`plot_tarp()`](https://pedroliman.github.io/neural.sbi/reference/plot_tarp.md)
  : Plot TARP expected coverage
- [`plot_posterior_predictive()`](https://pedroliman.github.io/neural.sbi/reference/plot_posterior_predictive.md)
  : Plot posterior predictive checks

## Benchmark tasks

Reference inference problems shared by the tests and benchmarks.

- [`task_gaussian_linear()`](https://pedroliman.github.io/neural.sbi/reference/tasks.md)
  [`task_two_moons()`](https://pedroliman.github.io/neural.sbi/reference/tasks.md)
  [`task_slcp()`](https://pedroliman.github.io/neural.sbi/reference/tasks.md)
  [`task_sir()`](https://pedroliman.github.io/neural.sbi/reference/tasks.md)
  : Benchmark tasks
