# NPE head-to-head against Python `sbi`

Run September 2026 with `neuralsbi` 0.6.71 and `sbi` 0.27.0 (torch 2.14, CPU).
Protocol and scripts: `inst/benchmarks/`, files `01`-`04`. Raw numbers:
`comparison_gaussian_linear_mdn.csv`, `comparison_gaussian_linear_maf.csv`,
`comparison_two_moons_mdn.csv`, `comparison_two_moons_maf.csv` beside this
file.

## What was run

Two tasks, two estimators, 10,000 simulations each, matching the roadmap's
M3 scope.

`task_gaussian_linear()`: prior `N(0, 0.1 I)` in 10 dimensions, likelihood
`x | theta ~ N(theta, 0.1 I)`, conjugate so `reference_posterior()` is exact.
2 observations, fresh draws from the prior.

`task_two_moons()`: prior `Uniform([-1, 1]^2)`, simulator maps `theta` onto a
noisy crescent whose posterior is bimodal for a typical observation. No
closed form, so only `neuralsbi` and `sbi` are compared to each other. 5
observations.

Both estimators, `mdn` and `maf`, trained on their own defaults (`sbi`'s
`NPE(density_estimator=...)`, `neuralsbi`'s `npe(density_estimator = ...)`):
MAF 5 transforms/50 hidden units, MDN 10 components/two 50-unit hidden
layers, batch 200, lr 5e-4, 10% validation, early stopping after 20 epochs
without improvement. Both sampled 10,000 draws per observation directly --
NPE samples the posterior in a forward pass, no MCMC on either side. Neither
side was tuned. `01_generate_data.R` writes the simulations once; both
implementations train on the exact same `(theta, x)` rows, so neither sees a
draw the other did not.

`gaussian_linear` ran at 2 observations rather than the script's default of
5. `c2st()`'s MLP classifier costs `O(d)` in its hidden-layer width (`10 * d`
units per layer, matching `sbibm`), and at `d = 10` a single 10000-vs-10000
comparison measured 956 seconds on this container's CPU (`d = 5` measures 373
seconds in `docs/verification-roadmap.md`). Five observations times three
comparisons per observation (`ours` vs `sbi`, `ours` vs reference, `sbi` vs
reference) would have cost several hours; two observations still shows
whether the two implementations agree, at a fraction of the wall clock.
`two_moons`, at `d = 2`, ran the full default of 5 observations -- a
10000-vs-10000 comparison there costs under a minute.

## Results

### gaussian_linear (dim 10, analytic reference)

| obs | estimator | c2st ours vs sbi | c2st ours vs ref | c2st sbi vs ref | max mean diff | max sd diff |
|---|---|---|---|---|---|---|
| 1 | mdn | 0.582 | 0.558 | 0.552 | 0.079 | 0.015 |
| 2 | mdn | 0.573 | 0.540 | 0.556 | 0.068 | 0.015 |
| 1 | maf | 0.535 | 0.522 | 0.518 | 0.056 | 0.017 |
| 2 | maf | 0.527 | 0.509 | 0.521 | 0.051 | 0.025 |

Every number is under the 0.60 bar, for both estimators, against `sbi` and
against the exact posterior. `max_mean_diff` and `max_sd_diff` are absolute
differences in posterior mean/sd against an exact posterior sd of `sqrt(1 /
(1/0.1 + 1/0.1)) = 0.224` -- so MDN's worst mean error (0.079) is about a
third of a standard deviation, and MAF's (0.056) about a quarter. MAF reads
closer to the coin-flip (0.50) than MDN on every column, which tracks: `sbi`
defaults to MAF, and it is the denser-tested path in both packages.

### two_moons (dim 2, no closed form)

| obs | estimator | c2st ours vs sbi | max mean diff | max sd diff |
|---|---|---|---|---|
| 1 | mdn | 0.643 | 0.013 | 0.005 |
| 2 | mdn | 0.621 | 0.015 | 0.008 |
| 3 | mdn | 0.657 | 0.005 | 0.003 |
| 4 | mdn | 0.647 | 0.018 | 0.002 |
| 5 | mdn | 0.625 | 0.005 | 0.001 |
| 1 | maf | 0.597 | 0.016 | 0.012 |
| 2 | maf | 0.587 | 0.042 | 0.016 |
| 3 | maf | 0.618 | 0.002 | 0.009 |
| 4 | maf | 0.596 | 0.030 | 0.017 |
| 5 | maf | 0.657 | 0.004 | 0.014 |

MDN misses the 0.60 bar at all 5 observations (0.62-0.66). MAF misses it at
2 of 5 (0.618, 0.657) and clears it narrowly at the other 3 (0.587-0.597).
Neither failure shows up in the moments: every mean difference is under 0.04
and every sd difference under 0.02, against a prior that spans `[-1, 1]` in
each dimension. `c2st()`'s classifier is an MLP, not the linear probe the
package used before it was aligned with `sbibm` (see
`docs/verification-roadmap.md`'s "C2ST aligned with `sbibm`" note): it can
separate two sample sets that agree on both first moments if the *shape*
differs, and two-moons' crescent is exactly the kind of shape a coordinatewise
mean/sd comparison is blind to. That is what the moments-vs-C2ST split is
doing here -- not two posteriors that clearly disagree, but two posteriors
close enough on every simple summary that only the shape comparison catches
the difference.

## Reading the acceptance criterion

`inst/benchmarks/README.md` sets C2ST <= 0.60 against `sbi`, and against the
analytic reference where one exists. **`gaussian_linear` passes outright**,
both estimators, both against `sbi` and against the exact posterior --
this is the easy case the bar was written for: a well-behaved, unimodal
Gaussian posterior that a forward-pass sampler recovers cleanly. **`two_moons`
does not clear the bar.** MDN misses at every observation; MAF is on the
line, passing 3 of 5 and missing 2 by 0.02-0.06.

That split is itself informative. `two_moons`' posterior is bimodal and
crescent-shaped -- the standard stress test in `sbibm` for exactly this
reason. MDN's Gaussian-mixture parameterization has to cover a curved,
non-convex region with axis-free but still elliptical components; MAF's
autoregressive flow can warp space more freely and lands closer to `sbi`'s
own MAF on every observation. Neither implementation has a reference to be
scored against on this task, so "does not clear 0.60" does not mean either
one is wrong -- it means two independently initialized, independently
optimized fits of a hard multimodal target are distinguishable at the C2ST's
resolution, which two default MAF/MDN fits of a bimodal density are not
guaranteed to avoid. `inst/benchmarks/two_moons_calibration.R` already checks
`neuralsbi`'s own two-moons fit (an NSF, not MDN/MAF) against SBC/TARP rather
than against `sbi`, and finds it calibrated -- consistent with the two moons
failure here being about MDN/MAF's fit to a hard shape rather than a
correctness problem in either package.

## Cost

Training itself is fast on both sides: MAF/MDN converge in tens to a few
hundred epochs at these simulation counts (`sbi`'s two_moons MAF took 181
epochs, `neuralsbi`'s 212; both `gaussian_linear` fits converged in under 65
epochs), each a matter of seconds to low minutes on this container's CPU. The
wall-clock cost is almost entirely `c2st()`'s MLP classifier: measured
directly on two 10000-draw standard-normal sets at `d = 10` (the case this
run needed), one `c2st()` call took 956 seconds. That is why
`gaussian_linear` ran at 2 observations instead of the script's default of 5
-- see "What was run" above.

## Reproducing

```sh
cd inst/benchmarks
Rscript 01_generate_data.R --task gaussian_linear --n 10000 --n_obs 2 --seed 42
python  02_run_sbi_python.py --task gaussian_linear --estimator maf --n_samples 10000 --seed 42
Rscript 03_run_neuralsbi.R --task gaussian_linear --estimator maf --n_samples 10000 --seed 42
Rscript 04_compare.R --task gaussian_linear --estimator maf

Rscript 01_generate_data.R --task two_moons --n 10000 --n_obs 5 --seed 42
python  02_run_sbi_python.py --task two_moons --estimator maf --n_samples 10000 --seed 42
Rscript 03_run_neuralsbi.R --task two_moons --estimator maf --n_samples 10000 --seed 42
Rscript 04_compare.R --task two_moons --estimator maf
```

Swap `--estimator maf` for `--estimator mdn` for the other estimator. Python
environment: `pip install torch --index-url
https://download.pytorch.org/whl/cpu` first, then `pip install sbi`, so `sbi`
does not pull a CUDA build of torch onto a CPU-only machine.
