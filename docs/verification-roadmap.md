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
- [ ] Better leakage handling (truncated proposals; log-prob normalization tests).
- [x] `summary()` methods (`nsbi_npe`, `nsbi_posterior`, `nsbi_samples`) and
      `as.data.frame.nsbi_samples()` tidy accessor (`R/summaries.R`).
- [x] `plot_coverage()` — nominal vs empirical coverage with Monte-Carlo band
      (`R/plotting.R`). Still open: TARP-style coverage; posterior-predictive
      plots.
- [ ] Vignettes: an applied end-to-end case study (e.g. SIR epidemic).
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
- [ ] **NSF** (rational-quadratic spline autoregressive flow) behind
      `density_estimator = "nsf"`. Plan: reuse `made_masks`; MADE outputs
      `3K - 1` spline parameters per dimension (K bin widths, K heights,
      K - 1 interior derivatives, softmax/softplus-constrained, linear tails
      outside `[-B, B]`, B ≈ 3); same stack/permutation structure as MAF.
- [ ] Verify each against analytic + `sbi` on SLCP and Two Moons.

### v0.4 — Embedding networks & structured data

- Optional learned summary network (MLP/CNN/RNN) mapping raw `x` → features,
  trained jointly. Enables time series / image-like observations.
- Standardization moves to embedding output.

### v0.5 — Sequential / multi-round NPE

- **NPE-C** (atomic proposal correction) and/or **NPE-A** (analytic correction)
  for simulation-efficient inference targeting a specific `x_o`.
- Proposal bookkeeping across rounds; importance-corrected loss.
- Verify simulation efficiency vs. single-round on fixed budgets.

### v0.6+ — Breadth

- Other families: NLE, NRE (ratio estimation) behind the same API.
- Ensemble posteriors; misspecification diagnostics; restriction estimators.
- Performance: GPU via torch, batched simulators, parallel simulation.

---

## Part C — Milestone checklist

- [x] M0 Pilot: linear-Gaussian analytic parity (torch-free) + MDN parity.
- [ ] M1 CI green with cached libtorch; Level-1 suite enforced.
- [ ] M2 Two Moons bimodality + SBC calibration demonstrated.
- [ ] M3 `sbi` head-to-head harness; Gaussian Linear & Two Moons C2ST ≤ 0.60.
- [ ] M4 MAF/NSF estimators; SLCP parity with `sbi`.
- [ ] M5 Embedding nets; a structured-data case study.
- [ ] M6 Sequential NPE-C; efficiency parity with `sbi`.
- [ ] M7 CRAN-ready: full docs, vignettes, `R CMD check` clean.

---

## Part D — Risks & mitigations

| Risk | Mitigation |
|---|---|
| libtorch install friction for users | ship exact `linear_gaussian` fallback; clear install docs; skip tests gracefully |
| MDN mode collapse / unstable training | validation early-stopping, restarts, standardization, full-cov components |
| Leakage with tight bounded priors | rejection + renormalized `log_prob`; roadmap: truncated proposals |
| Divergence from `sbi` internals | match z-scoring, optimizer, defaults; document any intentional differences |
| Flow correctness (v0.3) | test invertibility & log-det numerically; analytic + `sbi` parity gates |
