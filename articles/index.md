# Articles

### Learning the package

- [Introduction to neural simulation-based inference in
  R](https://pedroliman.github.io/neuralsbi/articles/intro-to-sbi.md):

  The basics of amortized, likelihood-free Bayesian inference: prior,
  simulator, observation, and how neuralsbi turns them into a posterior,
  illustrated with a simulator whose likelihood has no closed form.

- [Getting started with
  neuralsbi](https://pedroliman.github.io/neuralsbi/articles/neuralsbi.md):

  Train a neural posterior estimator from a prior and a simulator to
  perform amortized, likelihood-free Bayesian inference in R, then check
  that the fit is trustworthy.

- [Choosing a density
  estimator](https://pedroliman.github.io/neuralsbi/articles/density-estimators.md):

  Compare the conditional density estimators behind neural posterior
  estimation in R (mixture density network, masked autoregressive flow,
  and neural spline flow) and pick one for your simulator.

- [Checking the
  posterior](https://pedroliman.github.io/neuralsbi/articles/diagnostics.md):

  Validate a neural posterior in R with simulation-based calibration
  (SBC), expected coverage, TARP, and posterior predictive checks before
  you trust its credible intervals.

- [neuralsbi and pomp: two routes to an SIR
  posterior](https://pedroliman.github.io/neuralsbi/articles/sir-epidemic.md):

  A side-by-side comparison on a stochastic SIR epidemic: infer the
  transmission and recovery rates with pomp’s particle-filter MCMC and
  with neuralsbi’s neural posterior estimation, and check whether the
  two posteriors agree.
