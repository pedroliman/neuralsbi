# Changelog

## neuralsbi 0.2.4.9000 (development)

- Vignettes now show real output. They are *precomputed*: each
  vignette’s evaluated source lives in `vignettes/<name>.Rmd.orig`, and
  `vignettes/precompute.R` bakes it into a static `vignettes/<name>.Rmd`
  (results, printed values, and figures inlined). CI and pkgdown
  re-render that static Markdown with no torch at build time, so the
  expensive neural training runs once, locally, instead of on every
  build. Re-run `Rscript vignettes/precompute.R` after editing any
  `.Rmd.orig`.
- Two-moons calibration study
  (`inst/benchmarks/two_moons_calibration.R`): SBC, expected coverage,
  and TARP for a two-moons NSF fit, with figures written to
  `docs/figures/` (roadmap milestone M2).

## neuralsbi 0.2.3.9000

- Package website built with pkgdown, deployed from CI to
  <https://pedroliman.github.io/neural.sbi/>.
- Four vignettes that build on each other: getting started, choosing a
  density estimator, checking the posterior, and the SIR case study
  (which now also demonstrates
  [`npe_sequential()`](https://pedroliman.github.io/neural.sbi/reference/npe_sequential.md)).
  Removed a truncated duplicate of the SIR vignette.
- README rewritten to the standard terse form; authorship recorded in
  `DESCRIPTION` (Pedro Nascimento de Lima, with ORCID).

## neuralsbi 0.2.2.9000

- New
  [`npe_sequential()`](https://pedroliman.github.io/neural.sbi/reference/npe_sequential.md):
  multi-round NPE targeting a single observation via truncated-prior
  proposals (TSNPE, Deistler et al. 2022). Each round truncates the
  prior to the highest-probability region of the current posterior and
  retrains on all accumulated simulations; the standard NPE loss stays
  valid, so no importance correction is needed. Returns an `nsbi_snpe`
  fit that works with
  [`posterior()`](https://pedroliman.github.io/neural.sbi/reference/posterior.md),
  [`sample()`](https://pedroliman.github.io/neural.sbi/reference/sample.md),
  and the diagnostics, but is only valid at the targeted `x_obs`.
  Verified against the analytic linear-Gaussian posterior.

## neuralsbi 0.2.1.9000

- New
  [`tarp()`](https://pedroliman.github.io/neural.sbi/reference/tarp.md)
  diagnostic and
  [`plot_tarp()`](https://pedroliman.github.io/neural.sbi/reference/plot_tarp.md)
  (Lemos et al. 2023): a *joint* expected-coverage test using random
  reference points, complementing the per-parameter
  [`sbc()`](https://pedroliman.github.io/neural.sbi/reference/sbc.md)
  ranks. Detects posteriors with calibrated marginals but wrong
  correlation structure.
- New
  [`plot_posterior_predictive()`](https://pedroliman.github.io/neural.sbi/reference/plot_posterior_predictive.md):
  marginal predictive histograms with the observation marked; returns
  the observation’s predictive quantiles.
- Leakage correction is now under test: with a bounded prior, the
  renormalized
  [`log_prob()`](https://pedroliman.github.io/neural.sbi/reference/log_prob.md)
  integrates to one over the support and returns `-Inf` outside it
  (`test-posterior-normalization.R`).
- Fixed CI. `R CMD check` failed on three counts: the
  [`npe()`](https://pedroliman.github.io/neural.sbi/reference/npe.md)
  example required libtorch (it now uses the torch-free
  `linear_gaussian` estimator and runs unconditionally), the
  hand-maintained `npe.Rd`/`fit_mdn.Rd` usage sections had drifted
  behind the code (missing `n_restarts`, `clip_grad_norm`,
  `n_transforms`, and the `"maf"`/`"nsf"` options), and `CLAUDE.md` was
  not in `.Rbuildignore`. The `test-torch` job also failed because torch
  0.17 refuses a `TORCH_HOME` that does not exist; the workflow now
  creates it first.

## neuralsbi 0.2.0.9000

- Shared training engine for all neural estimators
  (`train_conditional_de()`): best-of-n restarts, learning-rate decay on
  plateau, gradient clipping, per-epoch loss history.
- Masked Autoregressive Flow (`density_estimator = "maf"`) and Neural
  Spline Flow (`"nsf"`, autoregressive rational-quadratic splines) join
  the MDN and the closed-form `linear_gaussian` baseline.
- Benchmark tasks
  ([`task_gaussian_linear()`](https://pedroliman.github.io/neural.sbi/reference/tasks.md),
  [`task_two_moons()`](https://pedroliman.github.io/neural.sbi/reference/tasks.md),
  [`task_slcp()`](https://pedroliman.github.io/neural.sbi/reference/tasks.md),
  [`task_sir()`](https://pedroliman.github.io/neural.sbi/reference/tasks.md))
  shared between tests and the `inst/benchmarks/` head-to-head benchmark
  harness.
- [`summary()`](https://rdrr.io/r/base/summary.html) methods,
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) tidy
  accessor,
  [`plot_coverage()`](https://pedroliman.github.io/neural.sbi/reference/plot_coverage.md).
- SIR applied case-study vignette.
- CI: `R CMD check` plus a `test-torch` job with cached libtorch.

## neuralsbi 0.1.0

- First pilot release: priors, single-round amortized
  [`npe()`](https://pedroliman.github.io/neural.sbi/reference/npe.md),
  `linear_gaussian` and MDN estimators, posterior sampling with leakage
  correction, SBC, expected coverage, C2ST, posterior-predictive checks,
  [`pairplot()`](https://pedroliman.github.io/neural.sbi/reference/pairplot.md),
  [`plot_sbc()`](https://pedroliman.github.io/neural.sbi/reference/plot_sbc.md).
