# Conditional density estimators

A conditional density estimator learns \\q\_\phi(\theta \mid x)\\. In
`neuralsbi` every estimator is trained in *standardized* space and
exposes two generics:

## Details

- `de_log_prob(de, theta, x)` – log density of `theta` given `x`

- `de_sample(de, x, n)` – draw `n` parameter vectors given a single `x`

Two estimators ship today:

- `"mdn"` – a Mixture Density Network (neural network -\> Gaussian
  mixture), the workhorse, requires the `torch` back end.

- `"linear_gaussian"` – a closed-form conditional Gaussian baseline
  (least-squares mean, residual covariance). No neural network, no
  `torch`. It is exact for linear-Gaussian simulators and doubles as a
  fast baseline and a regression-test oracle.
