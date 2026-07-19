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

### v0.2 — Robustness & ergonomics

- Restart-best-of-`n` training; learning-rate scheduling; gradient clipping.
- Better leakage handling (truncated proposals; log-prob normalization tests).
- `summary()` methods; tidy (`data.frame`) accessors for samples/diagnostics.
- `plot_coverage()`; TARP-style coverage; posterior-predictive plots.
- Vignettes: an applied end-to-end case study (e.g. SIR epidemic).
- CI with cached libtorch; `sbibm`-parity benchmark harness in `inst/benchmarks/`.

### v0.3 — Normalizing-flow density estimators

- **MADE** masked MLP → **MAF** (stacked MADE + permutations) and **NSF**
  (rational-quadratic spline coupling). These are `sbi`'s defaults and handle
  non-Gaussian posteriors (SLCP, etc.).
- Same `de_log_prob`/`de_sample` contract; selectable via
  `density_estimator = "maf" | "nsf"`.
- Verify each against analytic + `sbi` on SLCP and Two Moons.

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
