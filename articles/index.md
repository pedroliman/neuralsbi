# Articles

### Learning the package

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

- [Case study: inferring epidemic parameters
  (SIR)](https://pedroliman.github.io/neuralsbi/articles/sir-epidemic.md):

  A worked simulation-based inference case study in R: recover the
  transmission and recovery rates of a stochastic SIR epidemic model
  from a noisy incidence curve, with no tractable likelihood.
