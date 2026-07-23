# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this package is

`neuralsbi` is a **native R** implementation of neural simulation-based
inference, focused on Neural Posterior Estimation (NPE). It mirrors the workflow
of the Python [`sbi`](https://github.com/sbi-dev/sbi) package but is **not a
wrapper** — neural density estimators run on the R
[`torch`](https://torch.mlverse.org/) package (libtorch, no Python). Target
users are applied researchers, not ML engineers: sensible defaults, built-in
posterior diagnostics.

## Commands

Development uses standard R-package tooling. From the package root:

```r
# Load the package for interactive work (preferred over library() during dev)
devtools::load_all()

# Run the full test suite
devtools::test()

# Run a single test file
testthat::test_file("tests/testthat/test-maf.R")

# Regenerate NAMESPACE and man/*.Rd from roxygen comments (see caveat below)
devtools::document()

# Full R CMD check (what CI runs)
devtools::check()
```

From the shell, the CI-equivalent path (runs tests **inside** the package
namespace, so internal functions are visible to tests):

```sh
R CMD INSTALL --no-docs . && (cd tests && Rscript testthat.R)
```

Neural (torch) tests skip themselves automatically when libtorch is absent
(`tests/testthat/helper-torch.R`), so the suite runs everywhere; only the
`test-torch` CI job installs libtorch and exercises MDN/MAF/NSF end to end.

## Architecture

The pipeline is `prior + simulator → simulate_for_sbi() → standardize → density
estimator → npe fit → posterior() → sample()/log_prob()`. The pieces:

**The density-estimator contract is the central abstraction.** Every estimator
trains in **standardized** (z-scored) space and implements exactly two S3
methods:

- `de_log_prob(de, theta, x)` — length-n vector of log q(theta | x)
- `de_sample(de, x, n)` — n×dim matrix of draws given a single x

Standardization and its change-of-variables Jacobian live in `R/posterior.R`,
*not* in the estimators, so estimators stay simple and interchangeable. Adding a
new estimator means implementing these two methods and wiring a `fit_*` into the
`switch()` in `fit_density_estimator()` (`R/npe.R`). The estimators:

- `R/density_estimator.R` — `linear_gaussian`: closed-form conditional Gaussian,
  torch-free, **exact** for linear-Gaussian models. It is the regression oracle
  for the whole pipeline and lets most tests run without libtorch.
- `R/mdn.R` — Mixture Density Network (MLP → full-covariance Gaussian mixture).
- `R/flows.R` — MAF (masked autoregressive flow) + the shared MADE masking
  machinery (`made_masks`, `masked_linear`).
- `R/nsf.R` — NSF (autoregressive rational-quadratic spline flow). Reuses
  `made_masks` from `R/flows.R`. Note: `sbi`'s NSF uses coupling layers; this
  one is autoregressive (documented in-file).

**All neural estimators share one training loop:** `train_conditional_de()` in
`R/train.R` (train/val split, Adam, minibatching, early stopping, LR decay on
plateau via torch's `lr_reduce_on_plateau`, gradient clipping, best-of-n
restarts). Each `fit_*` passes a `build_net` closure and a `log_prob_fn`. Put
robustness/training changes here, once, rather than per estimator.

**Bounded priors leak.** A density estimator over unconstrained space places
mass outside a bounded prior's support. `R/posterior.R` handles this by (a)
rejection-sampling out-of-support draws and (b) renormalizing `log_prob` by the
estimated acceptance probability, returning `-Inf` outside support — matching
`sbi`'s strategy.

**Benchmark tasks are shared between tests and benchmarks.** `R/tasks.R` defines
`task_gaussian_linear` (with an analytic reference posterior), `task_two_moons`,
`task_slcp`, and `task_sir` as exported `nsbi_task` objects. Tests use them for
analytic-parity checks; `inst/benchmarks/` uses the same tasks for the
head-to-head comparison against Python `sbi`.

**Verification is layered** (see `docs/verification-roadmap.md`): Level 1 =
analytic ground truth in CI (linear-Gaussian is exact; flows checked with looser
tolerance); Level 2 = calibration (SBC, expected coverage) where no closed form
exists; Level 3 = head-to-head C2ST against Python `sbi` via `inst/benchmarks/`
(scripted, not run in CI). `docs/implementation-plan.md` covers the method and
module map; `docs/verification-roadmap.md` Part E is the live handoff/next-steps
section — **keep it current as you work.**

## R package conventions

- **`sample()` is an S3 generic** here (`R/generics.R`), masking `base::sample`.
  Posterior draws come from `sample.nsbi_posterior`. Keep the generic contract
  intact when touching sampling.
- **NAMESPACE and `man/*.Rd` are currently hand-maintained** (the files say so
  at the top). Roxygen comments are the source of truth in intent; when you add
  an exported function, either run `devtools::document()` or add the matching
  `export()`/`S3method()` line and a hand-written `.Rd` consistent with the
  others. Every export must have a man page — verify before committing.
- Keep the hard dependency surface minimal: `torch`, `testthat`, `knitr`,
  `rmarkdown` are all **Suggests**, never Imports. Neural code must degrade
  gracefully (skip, or point to `linear_gaussian`) when torch is unavailable —
  guard with `require_torch()` / the test helper, and never define a torch
  object at package-load time (wrap `nn_module()` construction in a function).
- Use `roxygen2` markdown docstrings; match the existing voice (see below).

## Git & release workflow

- **Never commit straight to `main`.** Branch for each task, push the branch,
  and open a PR — even for small changes. Prefer opening a GitHub **issue**
  first and linking the PR to it.
- **Increment the patch version often** (`DESCRIPTION` `Version:`) and note the
  change; the dev version carries a `.9000` suffix (currently `0.2.0.9000`).
  **Tag released versions** (`vX.Y.Z`) — there are no tags yet, so establish the
  habit.
- Commit in small, self-contained increments with a clear imperative subject and
  a body explaining the *why*. Run `devtools::check()` (or at least
  `devtools::test()`) before pushing.
- Keep `NEWS.md`/roadmap updated alongside code so the next contributor can pick
  up from the docs alone. **`NEWS.md` and version bumps are for substantive
  package code changes** — new features, behavior changes, bug fixes in `R/`.
  Cosmetic/tooling changes (README wording, the hex logo, pkgdown/SEO
  metadata, CI config) don't get a `NEWS.md` entry or a version bump.

## Writing style (docs, roxygen, commits, PRs, vignettes, comments)

Default to the register of a well-written R package — think the tone of the
`tidyverse`/`torch` docs: precise, plain, and confident. Write like a careful
human engineer, not a language model.

**Do:**

- Lead with the point. Short declarative sentences. Active voice.
- Explain *why* a thing exists or *why* a choice was made, not just what it does.
- Use concrete nouns and real function names. Prefer examples over adjectives.
- Vary sentence length naturally. Let technical precision carry the prose.
- In commit bodies, describe the change and its motivation as you would to a
  reviewing colleague.

**Do not** (these are LLM tells — avoid them):

- No filler openers: "In today's fast-paced world", "It's worth noting that",
  "It is important to note", "As we all know".
- No hype adjectives or adverbs: "powerful", "seamless", "robust" (as a
  buzzword), "cutting-edge", "leverage", "delve", "elevate", "unleash",
  "meticulous", "comprehensive", "effortlessly".
- No hollow summaries or self-congratulation: "In conclusion", "Overall, this
  provides a robust solution", "This ensures a smooth experience".
- No rule-of-three padding ("fast, reliable, and scalable") when one accurate
  word will do.
- No emoji in code, docs, or commits. No exclamation points in prose.
- Don't hedge everything ("might possibly perhaps"); state what is true and be
  direct about what is uncertain.
- Don't over-explain the obvious or restate the code in prose. Comment the
  non-obvious *why*.
- Avoid title-case headings that read like marketing; use sentence case.
- Don't hand-wrap prose with manual line breaks. Write each paragraph as one continuous line (README/`.Rmd`, `NEWS.md`, the `DESCRIPTION` Description field, `cran-comments.md`, roxygen, commit bodies) and let the editor or renderer soft-wrap. Automatic wrapping by tooling (e.g. pandoc's `github_document` default) is fine — that is the renderer's choice, not yours.
