# neuralsbi 0.2.3.9000 (development)

* Package website built with pkgdown, deployed from CI to
  https://pedroliman.github.io/neural.sbi/.
* Four vignettes that build on each other: getting started, choosing a
  density estimator, checking the posterior, and the SIR case study (which
  now also demonstrates `npe_sequential()`). Removed a truncated duplicate
  of the SIR vignette.
* README rewritten to the standard terse form; authorship recorded in
  `DESCRIPTION` (Pedro Nascimento de Lima, with ORCID).

# neuralsbi 0.2.2.9000

* New `npe_sequential()`: multi-round NPE targeting a single observation via
  truncated-prior proposals (TSNPE, Deistler et al. 2022). Each round truncates
  the prior to the highest-probability region of the current posterior and
  retrains on all accumulated simulations; the standard NPE loss stays valid,
  so no importance correction is needed. Returns an `nsbi_snpe` fit that works
  with `posterior()`, `sample()`, and the diagnostics, but is only valid at the
  targeted `x_obs`. Verified against the analytic linear-Gaussian posterior.

# neuralsbi 0.2.1.9000

* New `tarp()` diagnostic and `plot_tarp()` (Lemos et al. 2023): a *joint*
  expected-coverage test using random reference points, complementing the
  per-parameter `sbc()` ranks. Detects posteriors with calibrated marginals
  but wrong correlation structure.
* New `plot_posterior_predictive()`: marginal predictive histograms with the
  observation marked; returns the observation's predictive quantiles.
* Leakage correction is now under test: with a bounded prior, the renormalized
  `log_prob()` integrates to one over the support and returns `-Inf` outside
  it (`test-posterior-normalization.R`).
* Fixed CI. `R CMD check` failed on three counts: the `npe()` example required
  libtorch (it now uses the torch-free `linear_gaussian` estimator and runs
  unconditionally), the hand-maintained `npe.Rd`/`fit_mdn.Rd` usage sections
  had drifted behind the code (missing `n_restarts`, `clip_grad_norm`,
  `n_transforms`, and the `"maf"`/`"nsf"` options), and `CLAUDE.md` was not in
  `.Rbuildignore`. The `test-torch` job also failed because torch 0.17 refuses
  a `TORCH_HOME` that does not exist; the workflow now creates it first.

# neuralsbi 0.2.0.9000

* Shared training engine for all neural estimators (`train_conditional_de()`):
  best-of-n restarts, learning-rate decay on plateau, gradient clipping,
  per-epoch loss history.
* Masked Autoregressive Flow (`density_estimator = "maf"`) and Neural Spline
  Flow (`"nsf"`, autoregressive rational-quadratic splines) join the MDN and
  the closed-form `linear_gaussian` baseline.
* Benchmark tasks (`task_gaussian_linear()`, `task_two_moons()`,
  `task_slcp()`, `task_sir()`) shared between tests and the
  `inst/benchmarks/` head-to-head benchmark harness.
* `summary()` methods, `as.data.frame()` tidy accessor, `plot_coverage()`.
* SIR applied case-study vignette.
* CI: `R CMD check` plus a `test-torch` job with cached libtorch.

# neuralsbi 0.1.0

* First pilot release: priors, single-round amortized `npe()`, `linear_gaussian`
  and MDN estimators, posterior sampling with leakage correction, SBC, expected
  coverage, C2ST, posterior-predictive checks, `pairplot()`, `plot_sbc()`.
