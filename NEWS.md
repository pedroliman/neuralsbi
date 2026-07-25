# neuralsbi 0.4.1

* **Breaking: the simulator is now called once per parameter set and returns one simulated observation.** Most models a researcher already has work that way -- an ODE solve, an agent-based model, a call out to `pomp` or `deSolve` -- so meeting the old contract meant wrapping the model in an apply loop and thinking about column-major indexing before anything ran. Two signatures are accepted, decided from `formals(simulator)`: when every parameter name appears among the formals the parameters arrive by name, one scalar each (`function(mu, sigma) ...`); otherwise the whole named parameter vector goes to the first argument (`function(theta) ...`). A simulator returns a numeric vector of length `d`, a scalar, or a one-row matrix or data frame, and names on that output become the outcome names. See `?nsbi_simulator` for the contract and the migration examples.
* This also fixes a correctness bug introduced in 0.4.0. Chunking meant a simulator written for a single parameter set could return a vector whose length happened to match the chunk's row count, and that vector was accepted as one column of outcomes for several draws. `simulate_for_sbi(function(theta) c(theta[1] * 2, theta[2] * 2), prior, 100)` returned a 100 x 1 matrix of nonsense with no error and no warning; it now returns the 100 x 2 matrix it should.
* `npe()`, `simulate_for_sbi()`, `npe_sequential()`, `sbc()`, `tarp()` and `posterior_predictive()` gain `sim_args`, a named list forwarded to every simulator call. Observed data, a time grid, a fixed population size, a design matrix or a solver tolerance travel from the call site instead of being captured in a closure -- which also keeps them out of what gets serialized to a `future` worker. A list rather than `...` because every one of `npe()`'s formals sits before `...`, so R's partial matching would silently capture `x`, `theta`, `n` or `seed`.
* Simulations whose output contains `NA`, `NaN` or an infinite value are dropped, together with their parameters, with one warning per run reporting the count and the rate. A single non-finite value would otherwise poison the training loss and surface much later as a `NaN` validation loss. The count is recorded on the fit and shown by `print()`. Nothing left is an error. In `sbc()` and `tarp()` a failed draw removes the whole trial; in `posterior_predictive()` it reduces the number of predictive draws. Pre-computed `theta`/`x` passed to `npe()` are checked the same way. Note that dropping conditions on the simulator having succeeded: when failure depends on the parameters, the fit targets the posterior *given success*, so a high drop rate is a modelling signal.
* Random-number streams are now per simulation rather than per chunk, so a given seed produces the same simulations whatever the `future` plan and whatever the worker count. The 0.4.0 guarantee held only at a fixed chunk size.
* **`chunk_size` is gone**, from `npe()`, `simulate_for_sbi()`, `npe_sequential()`, `sbc()`, `tarp()` and `posterior_predictive()`, along with `options(neuralsbi.chunks)`. Chunking existed in 0.4.0 to make results reproducible across backends: the split had to depend on `n` alone, which meant users had to know about it and could change their draws by changing it. Per-simulation RNG streams give that guarantee outright, so what is left is ordinary parallel scheduling. Batches are now sized from the worker count, cannot affect a result, and are not a setting. Running sequentially there are no batches at all, just a loop.
* NSF's `n_bins` and `tail_bound` are explicit `npe()` arguments; `npe()` no longer takes `...`.
* New `save_npe()` / `load_npe()`. A fit whose estimator is `"maf"`, `"mdn"` or `"nsf"` holds a torch module, which is an external pointer: `saveRDS()` writes the pointer, the file reloads without complaint, and the first call that touches the network fails with `external pointer is not valid`. `save_npe()` writes the weights with `torch::torch_save()` and everything else as ordinary R objects, into one `.rds`; `load_npe()` rebuilds the network from the recorded architecture and restores them. An overnight fit can now be reloaded the next morning, which is what amortization was for. `posterior()` and `print()` also detect a fit that came back from `readRDS()` and say so, instead of failing later with a torch error.

# neuralsbi 0.4.0

* The simulator can now run in parallel. Declare a `future` plan -- `library(future); plan(multisession)` -- and every function that calls a simulator (`npe()`, `simulate_for_sbi()`, `npe_sequential()`, `sbc()`, `tarp()`, `posterior_predictive()`) spreads the work across workers. There is no new argument to pass and no parallel variant to call: with no plan declared everything runs sequentially as before, and `neuralsbi` mentions the two lines above once per session (`options(neuralsbi.parallel_hint = FALSE)` to silence it). Each chunk of parameters draws from its own L'Ecuyer-CMRG stream, so a given `set.seed()` produces the same simulations sequentially and on any number of workers. See `?nsbi_parallel`.
* Long-running work now reports progress with an ETA -- simulation and neural training alike, one progress step per simulation and per training epoch. With `progressr` installed, `neuralsbi` emits standard progressr updates, so `progressr::handlers()` and `with_progress()` control reporting; without it, a built-in bar needing no extra packages does the job. The training bar targets the epoch at which early stopping would fire and revises that target as the validation loss improves. See `?nsbi_progress`.
* `npe()`, `simulate_for_sbi()`, `npe_sequential()`, `sbc()`, `tarp()`, and `posterior_predictive()` gain a `chunk_size` argument controlling how many parameter rows go to the simulator per call. The default splits a run into about 64 chunks; because the split depends only on the number of simulations, results do not change with the number of workers. A simulator whose output must be produced in one call can set `chunk_size` to the full simulation budget.
* `future` and `progressr` are `Suggests`, not dependencies; `parallel` (base R) moves into `Imports`.
# neuralsbi 0.3.7

* The test suite now skips its plotting tests when `ggplot2`/`GGally` are not installed, instead of failing. Both are `Suggests`, so `R CMD check` under `_R_CHECK_FORCE_SUGGESTS_=false` -- the configuration CRAN uses on a machine without them -- previously hit 5 errors from `require_ggplot2()`. The new `skip_if_no_ggplot2()`/`skip_if_no_ggally()` helpers mirror the `skip_if_no_torch()` contract already used for the neural tests, so the suite runs everywhere.
* `vignette("sir-time-varying-beta")` reworks the epidemic model so that its posterior predictive tracks the observed case peaks instead of overshooting them. The introduction day and seed size are now inferred per state rather than fixed at two infections on 2020-01-21 (which left the epidemic's phase to demographic noise: replicate simulations at one fixed parameter differed by a factor of several thousand at the peak bin), the ascertainment prior is pinned by spring-2020 seroprevalence (case counts alone identify only the product of ascertainment and incidence, and the unconstrained fit slides toward vast, barely-ascertained epidemics), and the regime-duration parameterization drops a coordinate that the simulator's rescaling had left unidentified. The article reports the effective reproduction number *R*~e~(t) = beta(t) S(t) / (N gamma) computed from the simulator's own susceptible trajectories, adds a prior-predictive check before training, and summarizes the posterior predictive by its median rather than its mean.

# neuralsbi 0.3.6

* Parameters and outcomes can now be named. Name `prior_uniform()`'s `low`/`high` or `prior_normal()`'s `mean` (e.g. `c(beta = 0, gamma = 0)`), or attach `colnames()` to a simulator's output, and those names carry through `npe()`, `sample()`, `map_estimate()`, `sbc()`, `expected_coverage()`, and `posterior_predictive()` without any extra arguments. `plot_sbc()`, `plot_coverage()`, `pairplot()`, and `plot_posterior_predictive()` use the names for titles, legends, and facet labels; a name that happens to be valid R syntax (`"beta[1]"`, `"rho"`, `"sigma^2"`) renders as its plotmath symbol (Greek letters, sub/superscripts) instead of literal text.

# neuralsbi 0.3.5

* The SIR case study becomes a head-to-head comparison with the `pomp` package, `vignette("sir-epidemic")`. Both methods fit the same stochastic SIR epidemic: `pomp` via particle-filter MCMC (`pmcmc`), `neuralsbi` via neural posterior estimation. The vignette contrasts what each needs from the model — `pomp` a measurement density, `neuralsbi` only a simulator — overlays the two posteriors, scores their agreement with a C2ST, and confirms the neural fit with SBC. The comparison is precomputed, so `pomp` is needed only to regenerate the article, not to build or check the package.

# neuralsbi 0.3.4

* Plotting is now built on `ggplot2` and `GGally::ggpairs()` instead of base graphics. `pairplot()`, `plot_sbc()`, `plot_coverage()`, `plot_tarp()`, and `plot_posterior_predictive()` keep their signatures (`pairplot()` gains an `alpha` argument) but now build and print a `ggplot`/`ggmatrix` object, returned invisibly for further customization. `ggplot2` and `GGally` move to `Suggests`, following the same graceful-degradation pattern as `torch`: call any plotting function without them installed and you get an informative error, not a crash.
* New vignette, `vignette("intro-to-sbi")`: a short beginner tutorial covering the three ingredients (prior, simulator, observation), amortized training, and a first calibration check, using a g-and-k distribution simulator whose likelihood has no closed form.

# neuralsbi 0.3.2

* README is now generated from `README.Rmd`, so the usage example runs at render time and its output cannot drift from the code. The example is a plain linear regression: the posterior recovers the ground-truth coefficients, which the reader can check against ordinary least squares.

# neuralsbi 0.3.1

* CRAN resubmission fixes. Routed the `summary()` and `as.data.frame()`
  methods into the single `summaries` help topic (they had drifted into
  separate `.Rd` files with duplicated `\alias` entries, which also produced
  duplicate HTML anchors). Wrapped `theta_{<d}` in `\eqn{}` in the
  `made_masks` docs so the Rd no longer drops braces. Dropped the bare "NPE"
  acronym from the `DESCRIPTION` to avoid the spurious misspelling note.

# neuralsbi 0.3.0

* Defaults now match Python `sbi`, so a workflow reads the same in both
  packages and results can be cross-checked. Changes to `npe()` defaults:
  the density estimator is now `"maf"` (was `"mdn"`); MDN mixture components
  default to 10 (was 5); NSF spline bins default to 10 (was 8); the training
  batch size is 200 (was 100). `max_epochs` is raised to 2000 as a guard cap
  that early stopping (`patience = 20`) normally reaches first, mirroring
  `sbi`'s effectively-unbounded epoch budget. `lr`, `validation_fraction`,
  `patience`, `clip_grad_norm`, `n_transforms`, and `hidden` already matched.
  Pass any of these explicitly to recover the previous behavior.
* First CRAN submission. Dropped the development `.9000` version suffix,
  removed the redundant `Author`/`Maintainer` fields (now derived from
  `Authors@R`), and tidied the package title.

* Embedding networks (roadmap v0.4). `embedding_mlp()` builds a learned
  summary network that maps raw observations to a low-dimensional feature
  vector; pass it to `npe(..., embedding_net = )` and the MDN, MAF, and NSF
  estimators condition on the features instead of the raw data, training the
  embedding jointly. The estimators still take raw `x` at the `de_*` boundary
  (`dim_x` is unchanged), so sampling and `log_prob` route through the
  embedding automatically. Ignored, with a warning, by `linear_gaussian`.

# neuralsbi 0.2.4.9000 (development)

* Vignettes now show real output. They are *precomputed*: each vignette's
  evaluated source lives in `vignettes/<name>.Rmd.orig`, and
  `vignettes/precompute.R` bakes it into a static `vignettes/<name>.Rmd`
  (results, printed values, and figures inlined). CI and pkgdown re-render
  that static Markdown with no torch at build time, so the expensive neural
  training runs once, locally, instead of on every build. Re-run
  `Rscript vignettes/precompute.R` after editing any `.Rmd.orig`.
* Two-moons calibration study (`inst/benchmarks/two_moons_calibration.R`):
  SBC, expected coverage, and TARP for a two-moons NSF fit, with figures
  written to `docs/figures/` (roadmap milestone M2).

# neuralsbi 0.2.3.9000

* Package website built with pkgdown, deployed from CI to
  https://pedroliman.github.io/neuralsbi/.
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
