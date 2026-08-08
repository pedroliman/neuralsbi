# Conditional density estimators

A conditional density estimator learns \\q\_\phi(\theta \mid x)\\. In
`neuralsbi` every estimator is trained in *standardized* space and
exposes two generics:

## Details

- `de_log_prob(de, theta, x)` – log density of `theta` given `x`

- `de_sample(de, x, n)` – draw `n` parameter vectors given a single `x`

The contract is really `q(target | condition)`: it makes no assumption
about which of the two arguments is the parameter.
[`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) calls it
with `theta` as the target and `x` as the condition, learning the
posterior. [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md)
swaps the two, learning the likelihood instead, with the same estimators
and the same two generics.

Four estimators ship today:

- `"maf"` – a Masked Autoregressive Flow (Papamakarios et al., 2017), a
  stack of invertible autoregressive transforms with an exact
  change-of-variables density. This is the default, matching Python
  `sbi`, and requires the `torch` back end.

- `"nsf"` – a Neural Spline Flow (Durkan et al., 2019): the same
  autoregressive structure as the MAF, but with a monotonic
  rational-quadratic spline transform in place of MAF's affine one,
  which handles sharply non-Gaussian posteriors better. Requires
  `torch`.

- `"mdn"` – a Mixture Density Network (neural network -\> Gaussian
  mixture). Requires `torch`.

- `"linear_gaussian"` – a closed-form conditional Gaussian baseline
  (least-squares mean, residual covariance). No neural network, no
  `torch`. It is exact for linear-Gaussian simulators and doubles as a
  fast baseline and a regression-test oracle.

## References

Papamakarios, G., Pavlakou, T. and Murray, I. (2017). Masked
Autoregressive Flow for Density Estimation. *NeurIPS*.
[doi:10.48550/arXiv.1705.07057](https://doi.org/10.48550/arXiv.1705.07057)

Durkan, C., Bekasov, A., Murray, I. and Papamakarios, G. (2019). Neural
Spline Flows. *NeurIPS*.
[doi:10.48550/arXiv.1906.04032](https://doi.org/10.48550/arXiv.1906.04032)
