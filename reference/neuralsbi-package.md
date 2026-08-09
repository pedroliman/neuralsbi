# neuralsbi: Neural Simulation-Based Inference

A native R implementation of neural simulation-based inference (SBI)
methods. Given a prior over parameters and a simulator, 'neuralsbi'
trains a conditional neural density estimator for likelihood-free
Bayesian inference: Neural Posterior Estimation approximates the
posterior directly and samples it in a forward pass, Neural Likelihood
Estimation learns a surrogate likelihood that is sampled with MCMC,
handles repeated independent observations, and can be exported as 'Stan'
code for use inside a larger model, and Neural Ratio Estimation trains a
classifier for the likelihood-to-evidence ratio instead of either
density. Neural estimators run on the 'torch' back end. This package is
developed for applied researchers who want an R-native approachable SBI
interface with sensible defaults and built-in posterior diagnostics.

## See also

Useful links:

- <https://neuralsbi.pedrodelima.com/>

- <https://github.com/pedroliman/neuralsbi>

- Report bugs at <https://github.com/pedroliman/neuralsbi/issues>

## Author

**Maintainer**: Pedro Nascimento de Lima <plima@rand.org>
([ORCID](https://orcid.org/0000-0001-9057-198X))

Authors:

- Pedro Nascimento de Lima <plima@rand.org>
  ([ORCID](https://orcid.org/0000-0001-9057-198X))
