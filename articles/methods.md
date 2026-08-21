# Guide to neural SBI methods

This page maps each method in `neuralsbi` to the paper it comes from and
the function that runs it. For a worked example, read
[`vignette("neuralsbi")`](https://neuralsbi.pedrodelima.com/articles/neuralsbi.md)
instead.

## Key papers

- Cranmer et al. ([2020](#ref-cranmer2020frontier)) survey the field and
  set the NPE/NLE/NRE naming used here.
- Lueckmann et al. ([2021](#ref-lueckmann2021benchmarking)) compare
  these methods on a shared set of tasks.
- Deistler et al. ([2025](#ref-deistler2025practical)) give the applied
  workflow, from simulator and prior through training to diagnostics.

## Neural posterior estimation (NPE)

- Papamakarios and Murray ([2016](#ref-papamakarios2016npe)) introduced
  neural posterior estimation, training a conditional density estimator
  on prior draws so the fitted network is the posterior.
- [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md)
  implements the single-round amortized case. Conditioning on a new
  observation costs one forward pass through
  [`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md),
  then
  [`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md).

## Sequential neural posterior estimation

- Greenberg et al. ([2019](#ref-greenberg2019apt)) introduced automatic
  posterior transformation, which reweights the loss when draws come
  from a proposal rather than the prior.
- Deistler et al. ([2022](#ref-deistler2022tsnpe)) introduced truncated
  proposals, which restrict the prior to the current posterior’s
  highest-probability region and leave the plain NPE loss valid.
- [`npe_sequential()`](https://neuralsbi.pedrodelima.com/reference/npe_sequential.md)
  implements truncated proposals. The result is not amortized and holds
  only near the observation it was trained on.

## Neural likelihood estimation (NLE)

- Papamakarios et al. ([2018](#ref-papamakarios2019snl)) introduced
  sequential neural likelihood, learning the surrogate likelihood q(x
  \mid \theta) and recovering the posterior through Bayes’ rule.
- [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md)
  implements the single-round case, so repeated independent observations
  are a sum over one trained density. Evaluate it with
  [`log_lik()`](https://neuralsbi.pedrodelima.com/reference/log_lik.md)
  and sample it through
  [`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md).
- Neal ([2003](#ref-neal2003slice)) introduced slice sampling, which
  [`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md)
  uses by default.
  [`stan_code()`](https://neuralsbi.pedrodelima.com/reference/stan_export.md)
  writes the estimator out as Stan code instead.

## Neural likelihood-ratio estimation (NRE)

- Hermans et al. ([2019](#ref-hermans2020nre)) introduced amortized
  approximate likelihood-ratio estimation, training a binary classifier
  for r(\theta, x) = p(x \mid \theta) / p(x) rather than either density.
- Durkan et al. ([2020](#ref-durkan2020contrastive)) introduced the
  atomic contrastive objective, which scores the true parameter against
  contrasting draws from the same minibatch.
- [`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md)
  implements the atomic objective with 10 atoms per simulation by
  default. Evaluate it with
  [`log_ratio()`](https://neuralsbi.pedrodelima.com/reference/log_ratio.md)
  and sample it through
  [`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md).

## Density estimators

Select one through `density_estimator=` in
[`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) or
[`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md). A
closed-form `"linear_gaussian"` option is also available: it needs no
`torch`, is exact for linear-Gaussian models, and serves as a test
oracle.

### Mixture density network (MDN)

- Bishop introduced mixture density networks in the 1994 technical
  report [Mixture Density
  Networks](https://publications.aston.ac.uk/id/eprint/373/), which has
  no DOI.
- Papamakarios and Murray ([2016](#ref-papamakarios2016npe)) used an MDN
  as the density estimator for NPE.
- `npe(density_estimator = "mdn")` maps `x` to the weights, means, and
  full covariances of a Gaussian mixture over \theta. `n_components` and
  `hidden` size it.

### Masked autoregressive flow (MAF)

- Germain et al. ([2015](#ref-germain2015made)) introduced the MADE
  masking scheme for autoregressive density estimation.
- Papamakarios et al. ([2017](#ref-papamakarios2017maf)) introduced
  masked autoregressive flows, stacking affine autoregressive transforms
  to a standard-normal base.
- `npe(density_estimator = "maf")` is the default. `n_transforms` sets
  the depth.

### Neural spline flow (NSF)

- Durkan et al. ([2019](#ref-durkan2019nsf)) introduced neural spline
  flows, replacing the affine transform with a monotonic
  rational-quadratic spline for sharply non-Gaussian posteriors.
- `npe(density_estimator = "nsf")` implements the autoregressive form.
  `n_bins` and `tail_bound` control the spline.

## Diagnostics

### Simulation-based calibration (SBC)

- Talts et al. ([2018](#ref-talts2018sbc)) introduced simulation-based
  calibration, which ranks the true parameter among posterior draws over
  many prior draws. Calibrated posteriors give uniform ranks.
- [`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md) computes
  the ranks, one parameter at a time, and
  [`plot_sbc()`](https://neuralsbi.pedrodelima.com/reference/plot_sbc.md)
  displays them.

### Expected coverage

- Hermans et al. ([2021](#ref-hermans2022crisis)) showed that published
  SBI posteriors are often overconfident, and made coverage the check
  that exposes it.
- [`expected_coverage()`](https://neuralsbi.pedrodelima.com/reference/expected_coverage.md)
  turns the SBC ranks into nominal against empirical coverage at each
  credible level, and
  [`plot_coverage()`](https://neuralsbi.pedrodelima.com/reference/plot_coverage.md)
  displays them. Calibrated posteriors lie on the diagonal.

### TARP

- Lemos et al. ([2023](#ref-lemos2023tarp)) introduced tests of accuracy
  with random points, which measure coverage of distance-based credible
  regions around random reference points. This is a joint test, so it
  catches a posterior with calibrated marginals but wrong correlations.
- [`tarp()`](https://neuralsbi.pedrodelima.com/reference/tarp.md) runs
  the test and
  [`plot_tarp()`](https://neuralsbi.pedrodelima.com/reference/plot_tarp.md)
  displays it.

### Classifier two-sample test (C2ST)

- Lopez-Paz and Oquab ([2016](#ref-lopezpaz2017c2st)) introduced the
  classifier two-sample test, which trains a classifier to tell two sets
  of draws apart. Accuracy near 0.5 means they are indistinguishable.
- [`c2st()`](https://neuralsbi.pedrodelima.com/reference/c2st.md) uses
  cross-validated logistic regression on jointly z-scored draws. That
  classifier is linear, so it sees a shift in location but is close to
  blind to a difference in spread or dependence.

## References

Cranmer, Kyle, Johann Brehmer, and Gilles Louppe. 2020. “The Frontier of
Simulation-Based Inference.” *Proceedings of the National Academy of
Sciences* 117 (48): 30055–62. <https://doi.org/10.1073/pnas.1912789117>.

Deistler, Michael, Jan Boelts, Peter Steinbach, et al. 2025.
*Simulation-Based Inference: A Practical Guide*.
<https://doi.org/10.48550/arXiv.2508.12939>.

Deistler, Michael, Pedro J. Goncalves, and Jakob H. Macke. 2022.
*Truncated Proposals for Scalable and Hassle-Free Simulation-Based
Inference*. <https://doi.org/10.48550/arXiv.2210.04815>.

Durkan, Conor, Artur Bekasov, Iain Murray, and George Papamakarios.
2019. *Neural Spline Flows*.
<https://doi.org/10.48550/arXiv.1906.04032>.

Durkan, Conor, Iain Murray, and George Papamakarios. 2020. *On
Contrastive Learning for Likelihood-Free Inference*.
<https://doi.org/10.48550/arXiv.2002.03712>.

Germain, Mathieu, Karol Gregor, Iain Murray, and Hugo Larochelle. 2015.
*MADE: Masked Autoencoder for Distribution Estimation*.
<https://doi.org/10.48550/arXiv.1502.03509>.

Greenberg, David S., Marcel Nonnenmacher, and Jakob H. Macke. 2019.
*Automatic Posterior Transformation for Likelihood-Free Inference*.
<https://doi.org/10.48550/arXiv.1905.07488>.

Hermans, Joeri, Volodimir Begy, and Gilles Louppe. 2019.
*Likelihood-Free MCMC with Amortized Approximate Ratio Estimators*.
<https://doi.org/10.48550/arXiv.1903.04057>.

Hermans, Joeri, Arnaud Delaunoy, François Rozet, Antoine Wehenkel,
Volodimir Begy, and Gilles Louppe. 2021. *A Trust Crisis in
Simulation-Based Inference? Your Posterior Approximations Can Be
Unfaithful*. <https://doi.org/10.48550/arXiv.2110.06581>.

Lemos, Pablo, Adam Coogan, Yashar Hezaveh, and Laurence
Perreault-Levasseur. 2023. *Sampling-Based Accuracy Testing of Posterior
Estimators for General Inference*.
<https://doi.org/10.48550/arXiv.2302.03026>.

Lopez-Paz, David, and Maxime Oquab. 2016. *Revisiting Classifier
Two-Sample Tests*. <https://doi.org/10.48550/arXiv.1610.06545>.

Lueckmann, Jan-Matthis, Jan Boelts, David S. Greenberg, Pedro J.
Gonçalves, and Jakob H. Macke. 2021. *Benchmarking Simulation-Based
Inference*. <https://doi.org/10.48550/arXiv.2101.04653>.

Neal, Radford M. 2003. “Slice Sampling.” *The Annals of Statistics* 31
(3): 705–67. <https://doi.org/10.1214/aos/1056562461>.

Papamakarios, George, and Iain Murray. 2016. *Fast \epsilon-Free
Inference of Simulation Models with Bayesian Conditional Density
Estimation*. <https://doi.org/10.48550/arXiv.1605.06376>.

Papamakarios, George, Theo Pavlakou, and Iain Murray. 2017. *Masked
Autoregressive Flow for Density Estimation*.
<https://doi.org/10.48550/arXiv.1705.07057>.

Papamakarios, George, David C. Sterratt, and Iain Murray. 2018.
*Sequential Neural Likelihood: Fast Likelihood-Free Inference with
Autoregressive Flows*. <https://doi.org/10.48550/arXiv.1805.07226>.

Talts, Sean, Michael Betancourt, Daniel Simpson, Aki Vehtari, and Andrew
Gelman. 2018. *Validating Bayesian Inference Algorithms with
Simulation-Based Calibration*.
<https://doi.org/10.48550/arXiv.1804.06788>.
