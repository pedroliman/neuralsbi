# neuralsbi — Implementation Plan

Neural simulation-based inference in R, focused on **Neural Posterior
Estimation (NPE)**. A native R implementation modeled on the Python
[`sbi`](https://github.com/sbi-dev/sbi) package — *not* a wrapper around it.

---

## 1. Goals and scope

**Audience.** Applied researchers (epidemiology, ecology, economics,
psychology, engineering) who have a simulator and a prior and want a Bayesian
posterior without an intractable likelihood. They are comfortable in R and with
statistical modeling, but are *not* ML engineers.

**Design principles.**

1. **One obvious path.** `npe(prior, simulator)` → `posterior()` → `sample()`.
   Sensible defaults; every knob optional.
2. **Native R.** Neural density estimators are built on the
   [`torch`](https://torch.mlverse.org/) R package (libtorch C++ bindings, no
   Python). A pure-R closed-form estimator is always available as a fallback and
   oracle.
3. **Honest uncertainty.** Diagnostics (SBC, coverage, C2ST, posterior
   predictive) are first-class, not an afterthought.
4. **Verifiable.** Every method is checked against an analytic ground truth
   where one exists, and against Python `sbi` on shared benchmarks.

**In scope for v0.x:** amortized single-round NPE; MDN and (roadmap) normalizing
-flow density estimators; priors, simulators, posteriors, diagnostics, plots.

**Out of scope (for now):** NLE, NRE, sequential/multi-round methods, embedding
nets for high-dimensional/structured data, flow/score matching. All are on the
roadmap.

---

## 2. What `sbi` does, and what we mirror

The Python `sbi` workflow is:

```python
inference = NPE(prior)
inference.append_simulations(theta, x).train()
posterior = inference.build_posterior()
samples = posterior.sample((10_000,), x=x_obs)
```

Core ideas we reproduce:

| `sbi` concept | `neuralsbi` equivalent |
|---|---|
| `BoxUniform`, torch distributions as priors | `prior_uniform()`, `prior_normal()`, `prior_custom()` |
| `simulate_for_sbi` | `simulate_for_sbi()` |
| `NPE(prior).append_simulations().train()` | `npe(prior, simulator, ...)` |
| density estimators: MDN, MAF, NSF | `"mdn"` (v0), `"maf"`/`"nsf"` (roadmap) |
| z-scoring of theta/x | internal standardizers |
| `build_posterior()` + `.sample()` / `.log_prob()` | `posterior()` + `sample()` / `log_prob()` |
| leakage correction for bounded priors | rejection sampling + acceptance-renormalized `log_prob` |
| `sbc`, `run_sbc`, coverage, `c2st` | `sbc()`, `expected_coverage()`, `c2st()` |
| `pairplot` | `pairplot()`, `plot_sbc()` |

### The method, precisely

NPE trains a conditional density estimator \(q_\phi(\theta\mid x)\) by minimizing
the expected negative log-likelihood over simulations drawn from the prior:

\[
\mathcal{L}(\phi) = -\mathbb{E}_{\theta\sim p(\theta),\,x\sim p(x\mid\theta)}
\big[\log q_\phi(\theta \mid x)\big].
\]

When the proposal equals the prior (single round), the minimizer satisfies
\(q_\phi(\theta\mid x) \to p(\theta\mid x)\) — the true posterior — for *any*
\(x\) in the support of the marginal. This is why single-round NPE is already
correct and *amortized*: one training run yields a posterior for every possible
observation. Multi-round variants (NPE-A/C) only change the proposal to spend
simulations near a specific \(x_o\), and require a correction term; those are
roadmap items, not correctness prerequisites.

---

## 3. Architecture

```
prior ─┐
       ├─> simulate_for_sbi ─> (theta, x) ─> standardize ─> density estimator ─> nsbi_npe
sim  ──┘                                                          │
                                                                  ▼
                                          posterior(fit, x_obs) ─> nsbi_posterior
                                                                  │
                        sample() / log_prob() / map_estimate()  ◄─┤
                        sbc() / expected_coverage() / c2st()    ◄─┘
```

### Modules (`R/`)

| File | Responsibility |
|---|---|
| `prior.R` | `nsbi_prior` objects: sample, log_prob, support checks |
| `simulator.R`* | (folded into `npe.R`) run simulators over prior draws |
| `standardize.R` | z-scoring transforms + Jacobians |
| `density_estimator.R` | generic interface + closed-form `linear_gaussian` |
| `mdn.R` | torch Mixture Density Network (MLP → Gaussian mixture) |
| `npe.R` | `npe()` trainer, `simulate_for_sbi()` |
| `posterior.R` | `nsbi_posterior`: sample, log_prob, MAP, leakage handling |
| `diagnostics.R` | `sbc`, `expected_coverage`, `c2st`, `posterior_predictive` |
| `plotting.R` | `pairplot`, `plot_sbc` (base graphics) |
| `generics.R` | `sample()` S3 generic (sbi-style ergonomics) |
| `utils.R` | matrix coercion, torch availability |

### The density-estimator contract

Every estimator trains in **standardized space** and implements two S3 methods:

```r
de_log_prob(de, theta, x)  # length-n vector of log q(theta | x)
de_sample(de, x, n)        # n x dim matrix of draws given a single x
```

Standardization and the change-of-variables Jacobian live in `posterior.R`, so
estimators stay simple and interchangeable. Adding a normalizing flow later
means implementing exactly these two methods.

### The MDN (a neural estimator; the default is the MAF, matching `sbi`)

- MLP trunk (`hidden` widths, ReLU) maps `x` → shared features.
- Three linear heads output, per mixture component `k`:
  - mixture logit \(\ell_k\),
  - mean \(\mu_k \in \mathbb{R}^p\),
  - lower-triangular Cholesky factor \(L_k\) (diagonal via softplus > 0) so that
    \(\Sigma_k = L_k L_k^\top\) is a valid **full** covariance.
- Log-density via `logsumexp` over components; Cholesky solve for the quadratic
  form (numerically stable, no explicit inverse).
- Training: Adam, minibatches, validation split, early stopping on validation
  loss — matching `sbi`'s defaults (batch 100/50, lr 5e-4, 10% validation,
  patience ~20).

Full covariance (vs. diagonal) matters: the linear-Gaussian benchmark has a
correlated posterior, and we want to recover it exactly enough to pass C2ST.

---

## 4. Correctness & numerical choices

- **Standardization** is essential for training stability and is inverted
  (with Jacobian) for `log_prob`.
- **Leakage / bounded priors.** A density estimator over unconstrained space
  places some mass outside a bounded prior. We (a) reject out-of-support draws
  when sampling and (b) renormalize `log_prob` by the estimated acceptance
  probability, returning `-Inf` outside support — the same strategy as `sbi`.
- **Cholesky solves**, `logsumexp`, and softplus-parameterized scales keep the
  MDN loss finite and well-conditioned.
- **Reproducibility.** `seed` threads through simulation and torch.

---

## 5. Public API (v0.1)

```r
# priors
prior_uniform(low, high);  prior_normal(mean, sd);  prior_custom(...)
sample_prior(prior, n);    within_support(prior, theta)

# inference
simulate_for_sbi(simulator, prior, n)
fit <- npe(prior, simulator, n_simulations = 1000,
           density_estimator = c("mdn", "linear_gaussian"), ...)

# posterior
post <- posterior(fit, x_obs = ...)
sample(post, n);  log_prob(post, theta);  map_estimate(post)

# diagnostics & plots
sbc(fit, simulator);  expected_coverage(sbc_result);  c2st(a, b)
posterior_predictive(post, simulator)
pairplot(samples, truth = ...);  plot_sbc(sbc_result)
```

See `docs/verification-roadmap.md` for how each piece is validated against
analytic truth and against Python `sbi`.

*`simulator.R` responsibilities currently live in `npe.R`; split out if it grows.
