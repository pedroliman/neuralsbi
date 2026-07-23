# Mixture Density Network (MDN) conditional density estimator

The MDN is one of the neural density estimators in `neuralsbi` (the
default is the MAF, matching Python `sbi`). A multilayer perceptron maps
the data `x` to the parameters of a Gaussian mixture over the parameters
\\\theta\\: mixture logits, component means, and (full) lower-triangular
Cholesky factors of each component covariance. Training minimizes the
negative log-likelihood of \\\theta\\ under the mixture, which – when
simulations are drawn from the prior – yields a direct amortized
approximation of the posterior \\p(\theta \mid x)\\.

## Details

A native R/`torch` implementation of the multivariate-Gaussian mixture
density network (Bishop, 1994).
