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

## Acceptance criteria (roadmap M3)

On `gaussian_linear` and `two_moons` at 10k simulations:
C2ST(neuralsbi, sbi) <= 0.60, and both within C2ST <= 0.60 of the
reference posterior where one exists.

## File formats

- `data/<task>/theta.csv`, `data/<task>/x.csv` — one row per simulation, no header.
- `data/<task>/x_obs.csv` — one row per observation.
- `results/<task>/<impl>_<estimator>_obs<i>.csv` — posterior draws, one row per draw.

Python environment: `pip install sbi pandas` (sbi >= 0.22).
