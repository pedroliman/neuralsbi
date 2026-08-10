# Head-to-head benchmarks against Python `sbi`

Level-3 verification (see `docs/verification-roadmap.md`): train `neuralsbi`
and Python `sbi` on **identical simulations** and compare posteriors. Not run
in CI — run manually and commit the resulting metrics to `docs/benchmarks/`.

## Protocol

1. **Generate shared data** (R): draws `(theta, x)` from a task's prior and
   simulator plus a set of observations, written as CSVs to `data/<task>/`.

   ```sh
   Rscript 01_generate_data.R --task gaussian_linear --n 10000 --seed 42
   ```

2. **Train Python `sbi`** on those exact simulations; save posterior samples
   for each observation to `results/<task>/sbi_<estimator>_obs<i>.csv`.

   ```sh
   python 02_run_sbi_python.py --task gaussian_linear --estimator maf
   ```

3. **Train `neuralsbi`** on the same simulations; save samples to
   `results/<task>/neuralsbi_<estimator>_obs<i>.csv`.

   ```sh
   Rscript 03_run_neuralsbi.R --task gaussian_linear --estimator maf
   ```

4. **Compare** with C2ST, posterior mean/cov differences, and (where the task
   has an analytic reference) accuracy of both against ground truth:

   ```sh
   Rscript 04_compare.R --task gaussian_linear --estimator maf
   ```

## The NLE protocol (scripts 05–08)

`nle()` learns `q(x | theta)` from single observations and is then conditioned
on however many independent ones you have, so it needs its own data layout:
observation *sets* of several sizes drawn from one fixed parameter, rather than
one observation per row. Scripts `05`–`08` do that on `gaussian_linear`, whose
posterior stays conjugate for every set size, so both implementations are scored
against an exact answer rather than against each other.

```sh
Rscript 05_generate_data_nle.R --dim 5 --n 10000 --seed 42
python  06_run_sbi_nle.py      --estimator maf --n_samples 5000
Rscript 07_run_neuralsbi_nle.R --estimator maf --n_samples 5000
Rscript 08_compare_nle.R       --estimator maf
```

Both sides get the same simulations, the same estimator family and the same
sampler family (a vectorized slice sampler), and both run on their own
defaults. Tuning either one would answer a different question.

Read `c2st_*_vs_ref` next to the mean error and the sd ratio, not on its own:
`c2st()` here trains a logistic regression, which sees a shift in location and
is close to blind to a difference in spread (`?c2st` says so).

## The NRE protocol (scripts 09–12)

`nre()` learns the likelihood-to-evidence ratio `r(theta, x) = p(x | theta) /
p(x)` rather than a density, but it conditions on i.i.d. observation sets the
same way `nle()` does (see [`?nre`](../../man/nre.Rd)), so scripts `09`-`12`
reuse the NLE protocol's data layout on `task_gaussian_linear`, at `dim = 2`
by default to match the run this protocol reproduces (`docs/verification-roadmap.md`
Part E, #146).

```sh
Rscript 09_generate_data_nre.R --dim 2 --n 8000 --seed 42
python  10_run_sbi_nre.py      --classifier resnet
Rscript 11_run_neuralsbi_nre.R --classifier resnet
Rscript 12_compare_nre.R       --classifier resnet
```

`--classifier` takes `resnet`, `mlp` or `linear` on both sides (`sbi`'s three
`NRE`/`NRE_B` classifiers, matched by `nre()`'s architectures of the same
name); `neuralsbi`'s closed-form `logistic` classifier has no `sbi`
counterpart, so it only runs through `11`/`12`, scored against the reference
and, where a grid was generated, the analytic ratio.

`12_compare_nre.R` reports the same three numbers `08_compare_nle.R` does --
`c2st_ours_vs_ref`, `c2st_sbi_vs_ref`, `c2st_ours_vs_sbi` -- at every
observation-set size, which is the actual M3 bar for NRE and the piece #146
found missing: a `c2st()` between the two *posteriors*, not just a
correlation between the two *ratios*. When `09_generate_data_nre.R` also
wrote a grid (parameter dimension 2 or less, and 1 in `--n_obs`), `12` adds a
secondary comparison of the learned ratio itself against the analytic one --
correlation and centred RMSE, both series centred on their own mean since a
ratio is identified only up to an additive constant. That grid metric answers
"did the two classifiers learn the same shape"; the C2ST numbers answer "are
the two posteriors the same distribution", which is the one #146 asked for.

## Acceptance criteria (roadmap M3)

On `gaussian_linear` and `two_moons` at 10k simulations:
C2ST(neuralsbi, sbi) <= 0.60, and both within C2ST <= 0.60 of the
reference posterior where one exists. The NLE run holds itself to the same
0.60 against the conjugate reference, at every observation-set size. The NRE
run holds itself to the same 0.60, both `c2st_ours_vs_ref` and
`c2st_ours_vs_sbi`, at every observation-set size.

## File formats

- `data/<task>/theta.csv`, `data/<task>/x.csv` — one row per simulation, no header.
- `data/<task>/x_obs.csv` — one row per observation.
- `results/<task>/<impl>_<estimator>_obs<i>.csv` — posterior draws, one row per draw.

NLE runs use `data/gaussian_linear_iid/`: `x_obs_n<k>.csv` holds the first `k`
observations of one nested set, `theta_true.csv` the parameter they came from,
`meta.txt` the dimension followed by the set sizes. Draws land in
`results/gaussian_linear_iid/<impl>_<estimator>_n<k>.csv`.

NRE runs use `data/gaussian_linear_iid_nre/`, the same layout as the NLE
directory plus `grid_theta.csv` (only written when the parameter dimension is
2 or less and `1` is one of `--n_obs`): a fixed grid of parameter values
around the true one, conditioned on the first observation. `meta.txt` is the
dimension, the set sizes, and the grid size per dimension (`0` if no grid was
written). Draws land in
`results/gaussian_linear_iid_nre/<impl>_<classifier>_n<k>.csv`; grid log
ratios in `results/gaussian_linear_iid_nre/grid_log_ratio_<impl>_<classifier>.csv`.

Python environment: `pip install sbi pandas` (sbi >= 0.22; the NLE scripts were
run against 0.26.1). `sbi` depends on `nflows`, which does not build against
`setuptools >= 66` -- install `pip install "setuptools<66"` first if `pip
install sbi` fails while building `nflows`.
