# neuralsbi — Roadmap & Verification Plan

How we build out the package and **prove** it produces results in line with the
Python `sbi` package. Two intertwined tracks: **features** (what we implement)
and **verification** (how we know it's right).

---

## Part A — Verification strategy

We validate at three levels of increasing stringency.

### Level 1 — Analytic ground truth (unit/regression tests, run in CI)

For models with a closed-form posterior, we compare estimated posteriors to the
exact answer. No Python needed; fast; runs on every commit.

| Task | Ground truth | Metric & pass bar |
|---|---|---|
| **Linear Gaussian** \(x=\theta+\varepsilon\), \(\varepsilon\sim N(0,\sigma^2 I)\), Gaussian prior | Conjugate Gaussian posterior | posterior mean/sd within tol; C2ST < 0.55 vs analytic |
| **Gaussian, informative** (varying \(\sigma\), dims 1–5) | Conjugate Gaussian | C2ST < 0.6; coverage on diagonal |
| **Conjugate Beta–Bernoulli / Gamma–Poisson** | Analytic posterior | KS/C2ST vs analytic |

Status: **linear Gaussian implemented and passing** (see
`tests/testthat/test-linear-gaussian.R`, `test-mdn.R`). The `linear_gaussian`
estimator is *exact* here and acts as a torch-free oracle for the whole
pipeline; the MDN is checked against the same analytic target with looser tol.

### Level 2 — Calibration on models without closed-form posteriors

Correct posteriors are *calibrated*: the true parameter's rank within posterior
samples is uniform (SBC), and \(x\%\) credible intervals contain the truth \(x\%\)
of the time (expected coverage). These need no reference posterior at all.

| Task | Check | Pass bar |
|---|---|---|
| **Two moons** (bimodal) | SBC uniformity; visual bimodality | rank-hist uniform (p > 0.05); both modes recovered |
| **SIR / SIRD epidemic** (applied) | SBC + coverage | coverage within Monte-Carlo band |
| **Lotka–Volterra** | SBC + coverage | as above |

Tools already shipped: `sbc()`, `expected_coverage()`, `plot_sbc()`.

### Level 3 — Head-to-head against Python `sbi` (the headline claim)

On a shared, seeded benchmark suite, train both `neuralsbi` and Python `sbi` on
**identical simulations** and compare posteriors.

Protocol (scripted in `inst/benchmarks/`, not run in CI):

1. Fix a task, prior, seed. Generate `(theta, x)` once; save to disk.
2. Train Python `sbi` NPE (MDN, then MAF) on that exact dataset; save posterior
   samples for a fixed set of observations `x_o`.
3. Train `neuralsbi` on the same dataset; sample the same `x_o`.
4. Compare with:
   - **C2ST** between the two posteriors (target ≈ 0.5–0.55),
   - marginal & pairwise KS statistics,
   - posterior mean/cov difference,
   - **negative log-probability of held-out `theta`** under each posterior
     (proper scoring; should match within noise),
   - SBC rank histograms overlaid.
5. Report as a table + plots per (task × estimator × simulation budget).

**Reference tasks** (the `sbibm` benchmark set, which `sbi` itself uses): Gaussian
Linear, Gaussian Linear Uniform, Two Moons, SLCP, Bernoulli GLM, SIR, Lotka–
Volterra. We start with Gaussian Linear, Two Moons, SLCP, SIR.

**Acceptance criterion for "in line with `sbi`":** on Gaussian Linear and Two
Moons at ≥ 10k simulations, C2ST(neuralsbi, sbi) ≤ 0.60 and both within
C2ST ≤ 0.60 of the reference/analytic posterior, with overlapping SBC
histograms. Divergences are triaged (architecture, standardization, optimizer,
leakage handling) and documented.

### Continuous verification

- **CI (GitHub Actions):** `R CMD check` + Level-1 tests on every push; torch
  tests run on a job with libtorch cached, skipped elsewhere via
  `skip_if_no_torch()`.
- **Nightly/manual:** Level-2 calibration + Level-3 `sbi` comparison, artifacts
  (tables, plots) published to `docs/benchmarks/`.
- **Regression guard:** stored reference metrics; alert if C2ST/coverage drift.

---

## Part B — Feature roadmap

### v0.1 — Pilot (this release) ✅

- Priors: uniform, normal, custom.
- `npe()` amortized single-round NPE.
- Estimators: `linear_gaussian` (exact baseline/oracle), `mdn` (torch, full-cov
  Gaussian mixture).
- Posterior: sampling, `log_prob`, MAP, leakage correction.
- Diagnostics: SBC, expected coverage, C2ST, posterior predictive.
- Plots: `pairplot`, `plot_sbc`.
- Tests: prior, linear-Gaussian pipeline (analytic), MDN (analytic), diagnostics.
- Docs: implementation plan, this roadmap, README, getting-started vignette.

### v0.2 — Robustness & ergonomics (in progress)

- [x] Restart-best-of-`n` training; learning-rate decay on plateau; gradient
      clipping. All neural estimators share one training engine:
      `train_conditional_de()` in `R/train.R` (exposed via `npe(n_restarts =,
      clip_grad_norm =)`); per-epoch loss history stored in `fit$de$history`.
- [x] Better leakage handling: log-prob normalization tested
      (`test-posterior-normalization.R`: renormalized density integrates to 1
      over a bounded support, `-Inf` outside); truncated proposals shipped as
      TSNPE (`npe_sequential()`, see v0.5).
- [x] `summary()` methods (`nsbi_npe`, `nsbi_posterior`, `nsbi_samples`) and
      `as.data.frame.nsbi_samples()` tidy accessor (`R/summaries.R`).
- [x] `plot_coverage()` — nominal vs empirical coverage with Monte-Carlo band
      (`R/plotting.R`).
- [x] TARP joint coverage (`tarp()` in `R/diagnostics.R`, `plot_tarp()`;
      Lemos et al. 2023) and posterior-predictive plots
      (`plot_posterior_predictive()`). Tested with the linear-Gaussian oracle,
      including a miscalibration-detection case.
- [x] Vignettes: SIR applied case study (`vignettes/sir-case-study.Rmd`).
- [x] CI with cached libtorch: `.github/workflows/R-CMD-check.yaml` now has a
      `test-torch` job (installs libtorch, caches it, runs the full suite via
      `cd tests && Rscript testthat.R` so internals are visible to tests).
- [x] `sbibm`-parity benchmark harness in `inst/benchmarks/`: shared-data
      protocol scripted as `01_generate_data.R` → `02_run_sbi_python.py` →
      `03_run_neuralsbi.R` → `04_compare.R` (C2ST + moment diffs + analytic
      reference where available). Task definitions live in `R/tasks.R`
      (`task_gaussian_linear()`, `task_two_moons()`, `task_slcp()`, class
      `nsbi_task`) and are exported so tests and benchmarks share them.
      Still open: actually *running* the harness against Python `sbi` and
      committing the resulting metrics to `docs/benchmarks/` (needs a Python
      env with `pip install sbi`; see `inst/benchmarks/README.md`).

### v0.3 — Normalizing-flow density estimators (in progress)

- [x] **MADE** masked MLP → **MAF** (stacked MADE + order-reversal
      permutations, standard-normal base, identity-initialized transforms,
      clamped log-scales). Implemented in `R/flows.R`
      (`made_masks`/`made_module`/`maf_module`/`maf_forward`/`maf_inverse`),
      trained through the shared engine, selectable via
      `npe(density_estimator = "maf", n_transforms = )`.
      Tests: mask autoregressive invariants, forward/inverse round trip,
      identity-init log-prob equals standard normal, analytic linear-Gaussian
      parity (`tests/testthat/test-maf.R`).
- [x] **NSF** (rational-quadratic spline autoregressive flow) behind
      `npe(density_estimator = "nsf", n_transforms =, n_bins =, tail_bound =)`.
      Implemented in `R/nsf.R`: `rq_spline()` (batched monotonic RQ spline
      with linear tails outside `[-B, B]`, analytic inverse),
      `nsf_made_module()` (reuses `made_masks()`; emits `3K - 1` spline
      params/dim), `nsf_module()`/`nsf_log_prob_tensor()`/`nsf_inverse()`
      (same stack/reversal structure as MAF). Note: Python `sbi`'s NSF uses
      coupling layers; ours is autoregressive (documented in `R/nsf.R`).
      Tests: spline round trip + tail identity + log-det cancellation, full
      stack forward/inverse round trip, analytic linear-Gaussian parity
      (`tests/testthat/test-nsf.R`).
- [ ] Verify MAF/NSF against Python `sbi` on SLCP and Two Moons via
      `inst/benchmarks/` (harness is ready; needs a Python env with sbi).

### v0.4 — Embedding networks & structured data (in progress)

- [x] Optional learned summary network mapping raw `x` → features, trained
      jointly. `embedding_mlp(output_dim, hidden)` returns a torch-free
      `nsbi_embedding` spec; `npe(..., embedding_net = )` threads it into the
      MDN/MAF/NSF builders, which condition on the features. The embedding is
      a submodule of the fitted network (`embed_x()` applies it once per
      forward/inverse), so its parameters train with the estimator and travel
      in the `state_dict`; `de_log_prob`/`de_sample` still receive raw
      standardized `x` (`dim_x` unchanged). Implemented in `R/embedding.R`,
      wired through `R/mdn.R`, `R/flows.R`, `R/nsf.R`, `R/npe.R`. Tested in
      `tests/testthat/test-embedding.R` (spec validation torch-free; module
      shape, joint-training param registration, and linear-Gaussian parity
      with an embedded MAF under torch).
- [ ] CNN/RNN embeddings for image- and sequence-structured `x`; a
      structured-data case study (e.g. a time-series simulator).
- [ ] Standardize at the embedding output instead of the raw input (currently
      the embedding consumes the z-scored data; feature-space standardization
      would need a warm-up pass).

### v0.5 — Sequential / multi-round NPE (in progress)

- [x] **TSNPE** (truncated-prior proposals, Deistler et al. 2022):
      `npe_sequential(prior, simulator, x_obs, n_rounds, ...)` in
      `R/sequential.R`. Truncates the prior to the `1 - epsilon`
      highest-probability region of the current posterior via rejection,
      retrains on all accumulated rounds with the standard NPE loss (valid
      because every proposal is prior-proportional on its support). Tested
      against the analytic linear-Gaussian posterior, including bounded
      priors and per-round budgets (`test-sequential.R`).
- [ ] **NPE-C** (atomic proposal correction) and/or **NPE-A** (analytic
      correction); importance-corrected loss and proposal bookkeeping.
- [ ] Verify simulation efficiency vs. single-round on fixed budgets.

### v0.6+ — Breadth

- Other families: NLE, NRE (ratio estimation) behind the same API.
- Ensemble posteriors; misspecification diagnostics; restriction estimators.
- Performance: GPU via torch, batched simulators, parallel simulation.

---

## Part C — Milestone checklist

- [x] M0 Pilot: linear-Gaussian analytic parity (torch-free) + MDN parity.
- [~] M1 CI configured with cached libtorch (`test-torch` job) — needs one
      green run on GitHub to confirm.
- [~] M2 Two Moons bimodality test added (`test-two-moons.R`, torch-gated);
      SBC + coverage + TARP calibration study run and plotted
      (`inst/benchmarks/two_moons_calibration.R` → `docs/figures/`). Remaining:
      fold the figures into a short write-up / README section.
- [~] M3 `sbi` head-to-head harness scripted (`inst/benchmarks/`); running it
      and recording C2ST ≤ 0.60 still open.
- [~] M4 MAF and NSF estimators implemented + analytic parity tests; SLCP
      parity with `sbi` still open.
- [~] M5 MLP embedding net implemented (`embedding_mlp()`, `embedding_net`
      argument) + tests; CNN/RNN embeddings and a structured-data case study
      still open.
- [~] M6 Sequential NPE: TSNPE implemented + analytic parity test; NPE-C and
      efficiency parity with `sbi` still open.
- [ ] M7 CRAN-ready: full docs, vignettes, `R CMD check` clean.

---

## Part E — Handoff: current state & next actions

*Everything below is written so an agent (or human) with no other context can
pick up the work. Last updated for the 0.3.0 CRAN-prep pass (branch
`claude/sbi-cran-compliance-x1gtow`, July 2026): the `npe()` defaults were
aligned with Python `sbi` (density estimator `"maf"`, MDN `n_components = 10`,
NSF `n_bins = 10`, `batch_size = 200`, `max_epochs = 2000` as an early-stopping
guard), and the package was made CRAN-clean — `R CMD check --as-cran` reports
no package-level WARNINGs or NOTEs (only environmental ones: `torch`/qpdf/locale
absent, no network clock, badge 403 through the proxy). Version dropped its
`.9000` dev suffix to `0.3.0`; the redundant `Author`/`Maintainer` DESCRIPTION
fields now derive from `Authors@R`. The prior pass added MLP embedding networks
(`embedding_mlp()` + `embedding_net`), trained jointly inside MDN/MAF/NSF.*

### sbi default parity (0.3.0)

`npe()` and the `fit_*` estimators now default to the same hyperparameters as
Python `sbi`, so a workflow reads the same in both packages and results can be
cross-checked: estimator `"maf"`, `n_transforms = 5`, `hidden = c(50, 50)`,
MDN `n_components = 10`, NSF `n_bins = 10`/`tail_bound = 3`, `batch_size = 200`,
`lr = 5e-4`, `validation_fraction = 0.1`, `patience = 20`, `clip_grad_norm = 5`.
`max_epochs = 2000` is a guard cap; `sbi`'s epoch budget is effectively
unbounded and governed by early stopping, which the guard normally reaches
first. When changing a default, update the mirror in `fit_density_estimator()`
(the `%||%` fallbacks), the `fit_*` signatures, and the hand-written `.Rd`
`\usage` (codoc is checked).

### What exists right now

| Area | File(s) | State |
|---|---|---|
| Training engine | `R/train.R` | done; restarts, plateau LR decay, grad clipping, history |
| MDN | `R/mdn.R` | done, trains via shared engine |
| MAF | `R/flows.R` | done + tested (round trip, analytic parity) |
| NSF | `R/nsf.R` | done + tested; autoregressive (sbi uses coupling) |
| Tasks | `R/tasks.R` | gaussian_linear (analytic ref), two_moons, slcp, sir |
| Benchmarks vs sbi | `inst/benchmarks/01..04` | scripted, **never executed** |
| Summaries/tidy | `R/summaries.R` | done |
| Coverage plot | `plot_coverage()` in `R/plotting.R` | done |
| TARP coverage | `tarp()`, `plot_tarp()` | done + tested (calibrated & miscalibrated cases) |
| Posterior-predictive plot | `plot_posterior_predictive()` | done |
| Leakage normalization | tests in `test-posterior-normalization.R` | done |
| Sequential NPE (TSNPE) | `npe_sequential()` in `R/sequential.R` | done + analytic parity test; NPE-C open |
| Embedding net | `embedding_mlp()` in `R/embedding.R` | MLP done + tested; wired into MDN/MAF/NSF via `embedding_net`; CNN/RNN open |
| CI | `.github/workflows/R-CMD-check.yaml` | fixed (codoc drift, donttest example, TORCH_HOME); needs a green run on GitHub to confirm |
| NAMESPACE / man | hand-maintained | new exports have hand-written `.Rd`s |
| Website | `_pkgdown.yml`, `.github/workflows/pkgdown.yaml` | pkgdown site deployed to gh-pages; builds locally into `site/` (gitignored, `docs/` stays for these design docs); `pkgdown/strip-internal.R` removes CLAUDE.md from the output. GitHub Pages must be set to serve from the `gh-pages` branch once. |
| Vignettes | `vignettes/*.Rmd` (4) + `*.Rmd.orig` sources | neuralsbi (intro) → density-estimators → diagnostics → sir-epidemic. **Precomputed**: the evaluated source is `vignettes/<name>.Rmd.orig`; `vignettes/precompute.R` bakes it (with a working torch install) into a static `vignettes/<name>.Rmd` (results + figures inlined, figures under `vignettes/figures/`). CI and pkgdown re-render that static Markdown with no torch. Re-run `Rscript vignettes/precompute.R` after editing any `.Rmd.orig`; `.Rbuildignore` keeps the sources and the script out of the tarball. |
| Two-moons calibration | `inst/benchmarks/two_moons_calibration.R` | SBC + expected coverage + TARP on a two-moons NSF fit; figures in `docs/figures/two_moons_{sbc,coverage,tarp}.png` (M2) |

Key contract: every estimator implements `de_log_prob(de, theta, x)` and
`de_sample(de, x, n)` in **standardized** space (`R/density_estimator.R`);
`posterior.R` handles standardization, Jacobians, and leakage correction.
Neural estimators train via `train_conditional_de(build_net, log_prob_fn, ...)`.

### Environment notes (for a fresh container)

- Install R, then `install.packages(c("testthat", "torch"))` — set
  `USE_BUNDLED_LIBUV=1` if `fs` fails to compile — then
  `torch::install_torch()`.
- Run tests: `R CMD INSTALL --no-docs . && cd tests && Rscript testthat.R`
  (this runs tests inside the package namespace so internals are visible).
- torch-gated tests skip automatically without libtorch
  (`tests/testthat/helper-torch.R`).

### Next actions, in priority order

1. **CI is green on `main`** (M1 done) — the R-CMD-check workflow, including
   the `test-torch` job, last passed on the `main` merge commit. Longer term,
   consider generating NAMESPACE/man with roxygen2 so codoc drift can't recur.
2. **Run the sbi head-to-head** (finishes M3, headline claim). Needs
   `pip install sbi`. Follow `inst/benchmarks/README.md`: gaussian_linear
   and two_moons, estimators mdn + maf, 10k sims. Commit the comparison
   CSVs + a short summary to `docs/benchmarks/`. Acceptance: C2ST ≤ 0.60.
3. **Two-moons calibration study** (finishes M2): done via
   `inst/benchmarks/two_moons_calibration.R` — `sbc()` + `plot_coverage()` +
   `tarp()` on a two-moons NSF fit, figures saved to `docs/figures/`. TARP
   matters here: two-moons marginals hide the crescent structure that the
   joint test sees. Remaining: a short README/vignette section embedding the
   figures.
4. **SLCP with NSF** (finishes M4): `task_slcp()` exists; train NSF at 10k
   sims, compare to `sbi` via the harness. Expect this to stress leakage
   correction (uniform prior on [-3,3]^5).
5. **v0.2 leftovers**: none — log-prob normalization tests, TARP,
   posterior-predictive plots, the SIR vignette, and truncated proposals
   (via TSNPE) are all done.
6. **v0.4 embedding nets**: MLP done — `embedding_mlp()` + the `embedding_net`
   argument to `npe()` (`R/embedding.R`), trained jointly inside each
   estimator. Still open: CNN/RNN embeddings for image/sequence `x`, a
   structured-data case study, and (optionally) standardizing at the
   embedding output rather than the raw input.
7. **v0.5 sequential NPE**: TSNPE is done (`npe_sequential()`,
   `R/sequential.R`) and doubles as the truncated-proposal leakage
   handling. Still open: APT/NPE-C with the atomic loss correction, and a
   fixed-budget efficiency comparison (sequential vs. single-round) once a
   torch env is available.

### Known wrinkles / gotchas

- `torch_searchsorted(right = TRUE)` already returns the 1-based bin index
  for a (K+1)-edge grid (count of elements ≤ value) — do **not** add 1.
- NSF inversion: spline params must come from the partially reconstructed
  theta while the inverse acts on base-space values — see `nsf_apply()`'s
  `values` argument. A regression test covers this.
- The MDN mean-accuracy tests use tolerance 0.1–0.12; occasional seeds may
  be near the edge. If flaky in CI, bump sims, not tolerance.
- `sample()` masks `base::sample` (S3 generic in `R/generics.R`).
- DESCRIPTION `Version` is dev (`0.2.0.9000`); cut releases when M1–M3 are
  green.

---

## Part D — Risks & mitigations

| Risk | Mitigation |
|---|---|
| libtorch install friction for users | ship exact `linear_gaussian` fallback; clear install docs; skip tests gracefully |
| MDN mode collapse / unstable training | validation early-stopping, restarts, standardization, full-cov components |
| Leakage with tight bounded priors | rejection + renormalized `log_prob`; roadmap: truncated proposals |
| Divergence from `sbi` internals | match z-scoring, optimizer, defaults; document any intentional differences |
| Flow correctness (v0.3) | test invertibility & log-det numerically; analytic + `sbi` parity gates |
