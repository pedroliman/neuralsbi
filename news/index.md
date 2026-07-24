# Changelog

## neuralsbi 0.3.5

- The SIR case study becomes a head-to-head comparison with the `pomp`
  package,
  [`vignette("sir-epidemic")`](https://pedroliman.github.io/neuralsbi/articles/sir-epidemic.md).
  Both methods fit the same stochastic SIR epidemic: `pomp` via
  particle-filter MCMC (`pmcmc`), `neuralsbi` via neural posterior
  estimation. The vignette contrasts what each needs from the model —
  `pomp` a measurement density, `neuralsbi` only a simulator — overlays
  the two posteriors, scores their agreement with a C2ST, and confirms
  the neural fit with SBC. The comparison is precomputed, so `pomp` is
  needed only to regenerate the article, not to build or check the
  package.

## neuralsbi 0.3.4

- Plotting is now built on `ggplot2` and
  [`GGally::ggpairs()`](https://ggobi.github.io/ggally/reference/ggpairs.html)
  instead of base graphics.
  [`pairplot()`](https://pedroliman.github.io/neuralsbi/reference/pairplot.md),
  [`plot_sbc()`](https://pedroliman.github.io/neuralsbi/reference/plot_sbc.md),
  [`plot_coverage()`](https://pedroliman.github.io/neuralsbi/reference/plot_coverage.md),
  [`plot_tarp()`](https://pedroliman.github.io/neuralsbi/reference/plot_tarp.md),
  and
  [`plot_posterior_predictive()`](https://pedroliman.github.io/neuralsbi/reference/plot_posterior_predictive.md)
  keep their signatures
  ([`pairplot()`](https://pedroliman.github.io/neuralsbi/reference/pairplot.md)
  gains an `alpha` argument) but now build and print a
  `ggplot`/`ggmatrix` object, returned invisibly for further
  customization. `ggplot2` and `GGally` move to `Suggests`, following
  the same graceful-degradation pattern as `torch`: call any plotting
  function without them installed and you get an informative error, not
  a crash.
- New vignette,
  [`vignette("intro-to-sbi")`](https://pedroliman.github.io/neuralsbi/articles/intro-to-sbi.md):
  a short beginner tutorial covering the three ingredients (prior,
  simulator, observation), amortized training, and a first calibration
  check, using a g-and-k distribution simulator whose likelihood has no
  closed form.

## neuralsbi 0.3.2

- README is now generated from `README.Rmd`, so the usage example runs
  at render time and its output cannot drift from the code. The example
  is a plain linear regression: the posterior recovers the ground-truth
  coefficients, which the reader can check against ordinary least
  squares.

## neuralsbi 0.3.1

- CRAN resubmission fixes. Routed the
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods
  into the single `summaries` help topic (they had drifted into separate
  `.Rd` files with duplicated `\alias` entries, which also produced
  duplicate HTML anchors). Wrapped `theta_{<d}` in `\eqn{}` in the
  `made_masks` docs so the Rd no longer drops braces. Dropped the bare
  “NPE” acronym from the `DESCRIPTION` to avoid the spurious misspelling
  note.

## neuralsbi 0.3.0

- Defaults now match Python `sbi`, so a workflow reads the same in both
  packages and results can be cross-checked. Changes to
  [`npe()`](https://pedroliman.github.io/neuralsbi/reference/npe.md)
  defaults: the density estimator is now `"maf"` (was `"mdn"`); MDN
  mixture components default to 10 (was 5); NSF spline bins default to
  10 (was 8); the training batch size is 200 (was 100). `max_epochs` is
  raised to 2000 as a guard cap that early stopping (`patience = 20`)
  normally reaches first, mirroring `sbi`’s effectively-unbounded epoch
  budget. `lr`, `validation_fraction`, `patience`, `clip_grad_norm`,
  `n_transforms`, and `hidden` already matched. Pass any of these
  explicitly to recover the previous behavior.

- First CRAN submission. Dropped the development `.9000` version suffix,
  removed the redundant `Author`/`Maintainer` fields (now derived from
  `Authors@R`), and tidied the package title.

- Embedding networks (roadmap v0.4).
  [`embedding_mlp()`](https://pedroliman.github.io/neuralsbi/reference/embedding_mlp.md)
  builds a learned summary network that maps raw observations to a
  low-dimensional feature vector; pass it to
  `npe(..., embedding_net = )` and the MDN, MAF, and NSF estimators
  condition on the features instead of the raw data, training the
  embedding jointly. The estimators still take raw `x` at the `de_*`
  boundary (`dim_x` is unchanged), so sampling and `log_prob` route
  through the embedding automatically. Ignored, with a warning, by
  `linear_gaussian`.

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
  <https://pedroliman.github.io/neuralsbi/>.
- Four vignettes that build on each other: getting started, choosing a
  density estimator, checking the posterior, and the SIR case study
  (which now also demonstrates
  [`npe_sequential()`](https://pedroliman.github.io/neuralsbi/reference/npe_sequential.md)).
  Removed a truncated duplicate of the SIR vignette.
- README rewritten to the standard terse form; authorship recorded in
  `DESCRIPTION` (Pedro Nascimento de Lima, with ORCID).

## neuralsbi 0.2.2.9000

- New
  [`npe_sequential()`](https://pedroliman.github.io/neuralsbi/reference/npe_sequential.md):
  multi-round NPE targeting a single observation via truncated-prior
  proposals (TSNPE, Deistler et al. 2022). Each round truncates the
  prior to the highest-probability region of the current posterior and
  retrains on all accumulated simulations; the standard NPE loss stays
  valid, so no importance correction is needed. Returns an `nsbi_snpe`
  fit that works with
  [`posterior()`](https://pedroliman.github.io/neuralsbi/reference/posterior.md),
  [`sample()`](https://pedroliman.github.io/neuralsbi/reference/sample.md),
  and the diagnostics, but is only valid at the targeted `x_obs`.
  Verified against the analytic linear-Gaussian posterior.

## neuralsbi 0.2.1.9000

- New
  [`tarp()`](https://pedroliman.github.io/neuralsbi/reference/tarp.md)
  diagnostic and
  [`plot_tarp()`](https://pedroliman.github.io/neuralsbi/reference/plot_tarp.md)
  (Lemos et al. 2023): a *joint* expected-coverage test using random
  reference points, complementing the per-parameter
  [`sbc()`](https://pedroliman.github.io/neuralsbi/reference/sbc.md)
  ranks. Detects posteriors with calibrated marginals but wrong
  correlation structure.
- New
  [`plot_posterior_predictive()`](https://pedroliman.github.io/neuralsbi/reference/plot_posterior_predictive.md):
  marginal predictive histograms with the observation marked; returns
  the observation’s predictive quantiles.
- Leakage correction is now under test: with a bounded prior, the
  renormalized
  [`log_prob()`](https://pedroliman.github.io/neuralsbi/reference/log_prob.md)
  integrates to one over the support and returns `-Inf` outside it
  (`test-posterior-normalization.R`).
- Fixed CI. `R CMD check` failed on three counts: the
  [`npe()`](https://pedroliman.github.io/neuralsbi/reference/npe.md)
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
  ([`train_conditional_de()`](https://pedroliman.github.io/neuralsbi/reference/train_conditional_de.md)):
  best-of-n restarts, learning-rate decay on plateau, gradient clipping,
  per-epoch loss history.
- Masked Autoregressive Flow (`density_estimator = "maf"`) and Neural
  Spline Flow (`"nsf"`, autoregressive rational-quadratic splines) join
  the MDN and the closed-form `linear_gaussian` baseline.
- Benchmark tasks
  ([`task_gaussian_linear()`](https://pedroliman.github.io/neuralsbi/reference/tasks.md),
  [`task_two_moons()`](https://pedroliman.github.io/neuralsbi/reference/tasks.md),
  [`task_slcp()`](https://pedroliman.github.io/neuralsbi/reference/tasks.md),
  [`task_sir()`](https://pedroliman.github.io/neuralsbi/reference/tasks.md))
  shared between tests and the `inst/benchmarks/` head-to-head benchmark
  harness.
- [`summary()`](https://rdrr.io/r/base/summary.html) methods,
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) tidy
  accessor,
  [`plot_coverage()`](https://pedroliman.github.io/neuralsbi/reference/plot_coverage.md).
- SIR applied case-study vignette.
- CI: `R CMD check` plus a `test-torch` job with cached libtorch.

## neuralsbi 0.1.0

- First pilot release: priors, single-round amortized
  [`npe()`](https://pedroliman.github.io/neuralsbi/reference/npe.md),
  `linear_gaussian` and MDN estimators, posterior sampling with leakage
  correction, SBC, expected coverage, C2ST, posterior-predictive checks,
  [`pairplot()`](https://pedroliman.github.io/neuralsbi/reference/pairplot.md),
  [`plot_sbc()`](https://pedroliman.github.io/neuralsbi/reference/plot_sbc.md).
