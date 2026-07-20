# Dev Notes — Conversation Summary & Next Steps

_Last updated: 2026-07-20_

This file captures the state of the `neuralsbi` project, what was decided, and
what to do next. It is developer-facing scratch documentation, not user docs.

---

## 1. Goal

Build a **native R** package for **neural simulation-based inference (SBI)**,
focused on **Neural Posterior Estimation (NPE)**, modeled on the Python
[`sbi`](https://github.com/sbi-dev/sbi) package but **not a wrapper**. Neural
density estimators run on the R [`torch`](https://torch.mlverse.org/) package
(libtorch C++ bindings — no Python). Target users are **applied researchers**,
not ML engineers: easy defaults, built-in diagnostics.

---

## 2. What exists now (v0.1 pilot — shipped)

A working, tested pilot. All of the following is implemented, runs, and is
verified against analytic ground truth.

- **Priors**: `prior_uniform`, `prior_normal`, `prior_custom` — sampling,
  log-density, support checks (`R/prior.R`).
- **`npe()`** amortized single-round NPE trainer + `simulate_for_sbi()`
  (`R/npe.R`).
- **Density estimators** behind a common contract (`de_sample`, `de_log_prob`):
  - `"mdn"` — torch Mixture Density Network: MLP → **full-covariance** Gaussian
    mixture, logsumexp NLL, Adam + early stopping, z-scoring (`R/mdn.R`).
  - `"linear_gaussian"` — closed-form conditional Gaussian, torch-free, exact
    for linear-Gaussian models; baseline + regression oracle
    (`R/density_estimator.R`).
- **Posterior**: `posterior()`, `sample()`, `log_prob()`, `map_estimate()`;
  rejection + acceptance-renormalization for bounded priors (leakage handling)
  (`R/posterior.R`).
- **Diagnostics**: `sbc()`, `expected_coverage()`, `c2st()`,
  `posterior_predictive()` (`R/diagnostics.R`).
- **Plots**: `pairplot()`, `plot_sbc()` (base graphics, dependency-free)
  (`R/plotting.R`).
- **Docs**: `docs/implementation-plan.md`, `docs/verification-roadmap.md`,
  getting-started vignette, README with result figures.
- **Infra**: roxygen man pages, hand-written `NAMESPACE`, GitHub Actions
  `R-CMD-check` workflow.

### Verified (actually run in the dev container)
- **25 testthat tests pass** (torch MDN tests self-skip when libtorch absent).
- MDN **recovers the analytic linear-Gaussian posterior**, C2ST ≈ 0.53 vs.
  analytic (0.5 = indistinguishable) — `docs/figures/linear_gaussian_posterior.png`.
- MDN **recovers the bimodal two-moons posterior** —
  `docs/figures/two_moons_posterior.png`.
- `R CMD check` clean on code/docs/examples/tests. Remaining warnings are
  environmental only (container lacks `rmarkdown`/`pandoc` + UTF-8 locale; CI
  resolves them).

### Environment setup performed in the container
- Installed R 4.3.3 (`apt-get install r-base r-base-dev`).
- Installed R `torch` (compiled from source; libtorch downloaded) — **works**.
- Installed `testthat` (needed `apt libuv1-dev` for `fs`), `roxygen2` (needed
  `apt libxml2-dev` for `xml2`), and `S7`.

---

## 3. Decisions taken this conversation

1. **Package name** = `neuralsbi` (repo stays `neural.sbi`). Dotted package
   names cause `library()`/install friction. Revisit only if the user insists on
   matching the repo name.
2. **Backend** = R `torch` (native libtorch), never a Python bridge. Keep the
   torch-free `linear_gaussian` estimator as a permanent fallback + test oracle.
3. **NPE method** = single-round amortized NPE for v0.1 (theoretically complete:
   with prior as proposal, the trained `q(θ|x)` → true posterior). Sequential
   NPE-A/C are roadmap items, not correctness prerequisites.
4. **OO SYSTEM → migrate to S7** (user directive, this conversation). The pilot
   currently uses **S3**; refactor to **S7** (fall back to R6 only where mutable
   reference semantics are genuinely needed). See §5 for the plan.

---

## 4. `sbi` architecture reference (researched this conversation)

Grounding for the roadmap. Source tree of `sbi-dev/sbi` under `sbi/`:

| Module | Contents |
|---|---|
| `inference/` | Trainers (NPE/NLE/NRE) + base classes. |
| `neural_nets/` | `estimators/`, `embedding_nets/`, `net_builders/`, `factory.py`, `ratio_estimators.py`. |
| `posteriors/` | `DirectPosterior`, `MCMCPosterior`, `RejectionPosterior`, `VIPosterior`, `ImportanceSamplingPosterior`, `EnsemblePosterior`. |
| `samplers/` | MCMC / rejection / VI sampling utilities. |
| `diagnostics/` | `run_sbc`/`check_sbc`/`sbc_rank_plot`, `run_tarp`/`plot_tarp`, `c2st`. |
| `analysis/` | `pairplot` and result analysis. |
| `simulators/`, `utils/` | Simulator interfaces, helpers. |

**Factory model strings** (`neural_nets/factory.py`): mixture `mdn`; flows
`made`, `maf`, `maf_rqs`, `nsf`, `tabpfn`; mixed/categorical `mnle`, `mnpe`;
Zuko flows `zuko_{nice,maf,nsf,ncsf,sospf,naf,unaf,gf,bpf}`; score/flow-matching
via `posterior_score_nn` / `posterior_flow_nn`. Factories: `posterior_nn`
(NPE), `likelihood_nn` (NLE), `classifier_nn` (NRE), plus score/flow/marginal.
Key hyperparameters: `hidden_features`, `num_transforms`, `num_bins`,
`num_components`, `embedding_net`. **Z-scoring modes**: `independent`,
`structured`, `transform_to_unconstrained`, `none`.

**Algorithms implemented in `sbi`:**
- NPE: NPE(default, flow), **NPE_A** (MDN; Papamakarios & Murray 2016),
  **NPE_B** (Lueckmann 2017), **NPE_C/APT** (Greenberg 2019), TSNPE, **FMPE**
  (flow matching; Dax 2023), **NPSE** (score/diffusion; Geffner/Sharrock 2023–24).
- NLE: NLE (Papamakarios 2019), **MNLE** (mixed; Boelts 2022).
- NRE: NRE_A (Hermans 2020), NRE_B/SRE (Durkan 2020), NRE_C (Miller 2022),
  BNRE (Delaunoy 2022).
- Diagnostics: `run_sbc`, `check_sbc`, `run_tarp`, `c2st`; `RestrictionEstimator`.

Our v0.1 corresponds to **NPE_A** (MDN density estimator, single round).

---

## 5. NEXT STEPS

### 5.1 Immediate: refactor S3 → S7 (in progress, user-directed)

S7 (`library(S7)`, v0.2.2 installed) is the modern R OO system. Verified in the
container: `new_class`, properties, `validator`, inheritance, and
`new_generic`/`method<-` dispatch all work. **Gotcha:** `dim` is a **reserved
property name** — use `n_dim` / `param_dim` instead.

Proposed S7 class hierarchy (mirrors `sbi`):

```
Prior            (base)                properties: n_dim, lower, upper
  ├─ UniformPrior                      + (bounds)
  ├─ NormalPrior                       + mean, sd
  └─ CustomPrior                       + sample_fn, log_prob_fn
   generics: draw_prior(prior, n), prior_log_prob(prior, theta),
             in_support(prior, theta)

DensityEstimator (base)               properties: n_dim_theta, n_dim_x, fitted state
  ├─ MDN            (holds torch nn_module)
  └─ LinearGaussian (holds B, Sigma)
   generics: de_train(spec, theta, x), de_log_prob(de, theta, x),
             de_sample(de, x, n)

Trainer / NPE     (builder, mirrors sbi's NPE())   properties: prior, theta, x,
                                                    estimator, std_theta, std_x
   generics: append_simulations(npe, theta, x), train(npe, ...),
             build_posterior(npe, x_obs)

Posterior (base) → DirectPosterior    properties: estimator, prior, stds, default_x
   generics: sample(post, n, x), log_prob(post, theta, x), map(post, x)

Result value classes: PosteriorSamples, SBCResult, CoverageResult
```

Refactor tasks:
- [ ] `Depends: R (>= 4.1.0)` → add `Imports: S7`; `.onLoad` calls
      `S7::methods_register()`.
- [ ] Port `R/prior.R`, `R/density_estimator.R`, `R/mdn.R`, `R/npe.R`,
      `R/posterior.R` to S7 classes + generics. Keep torch nn_module as a
      property (reference held inside an S7 value object is fine).
- [ ] Keep the **functional convenience API** (`npe()`, `posterior()`,
      `sample()`, `log_prob()`, `sbc()`, `c2st()`, `pairplot()`) as thin
      wrappers so the user-facing surface and all 25 tests keep working.
- [ ] Also expose the **sbi-style class API**:
      `NPE(prior) |> append_simulations(θ,x) |> train() |> build_posterior()`.
- [ ] Decide `sample()`: keep as an S7 generic that masks `base::sample` with a
      `class_any` fallback delegating to `base::sample` (same tradeoff as the
      current S3 approach). Alternative: primary `draw()`/`sample_posterior()`.
- [ ] Regenerate roxygen (`@include` ordering matters for S7 class defs — add
      `Collate` or `@include` tags so classes are defined before methods).
- [ ] Re-run `R CMD check` + testthat; update vignette/README snippets.

### 5.2 Roadmap (feature + verification) — see `docs/verification-roadmap.md`

Prioritized milestones:
- **M1** CI green with cached libtorch; Level-1 analytic tests enforced.
- **M2** Two Moons bimodality + SBC calibration demonstrated (mostly done).
- **M3 (headline)** `sbi` head-to-head harness (`inst/benchmarks/`): train R
  `neuralsbi` and Python `sbi` NPE on **identical seeded simulations**
  (`sbibm` tasks: Gaussian Linear, Two Moons, SLCP, SIR); compare with C2ST,
  marginal KS, SBC overlays. Acceptance: C2ST(neuralsbi, sbi) ≤ 0.60 on Gaussian
  Linear + Two Moons at ≥10k sims.
- **M4** Normalizing-flow estimators: **MADE → MAF → NSF** (spline coupling),
  same `de_*` contract. This is `sbi`'s default family; needed for SLCP and
  heavy-tailed/non-Gaussian posteriors. Verify vs analytic + `sbi`.
- **M5** Embedding networks (MLP/CNN/permutation-invariant) for structured/
  high-dim `x`; z-scoring moves to embedding output.
- **M6** Sequential NPE-C (atomic proposal correction) / NPE-A analytic
  correction for simulation efficiency; add `RestrictionEstimator`.
- **M7** TARP diagnostic (`run_tarp`) alongside SBC; richer coverage plots.
- **M8** Breadth behind the same API: NLE, NRE; ensemble & VI posteriors;
  misspecification checks. CRAN-readiness.

### 5.3 Known rough edges to address
- Two-moons crescents are coarse at pilot settings (10k sims, 6 components,
  CPU). Cleaner modes need more sims/components or the flow estimators (M4).
- MDN training is single-threaded-ish on CPU; set `torch::torch_set_num_threads`
  and consider GPU support. `mdn_build_tril()` loops over lower-tri entries —
  fine for small `p`, revisit for large parameter dimension.
- Vignette needs `pandoc` to build; ensure CI image has it.

---

## 6. How to reproduce locally
```bash
Rscript -e 'library(testthat); library(neuralsbi); test_dir("tests/testthat")'
Rscript inst/examples/pilot_demo.R      # regenerates both figures
R CMD build neural.sbi && R CMD check neuralsbi_0.1.0.tar.gz
```
