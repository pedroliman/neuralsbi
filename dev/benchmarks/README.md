# Reproducing "Benchmarking Simulation-Based Inference" with neuralsbi

These scripts run the benchmark from Lueckmann, Boelts, Greenberg, Gonçalves and Macke (2021), *Benchmarking Simulation-Based Inference* (AISTATS), using `neuralsbi` in place of Python `sbi`, and tell you cell by cell whether we land where the paper landed.

The paper's grid is 10 tasks x 8 algorithms x 3 simulation budgets x 10 observations. We cover the two algorithms `neuralsbi` implements, NPE and NLE, across all 10 tasks and all three budgets (10³, 10⁴, 10⁵) and all 10 observations. That is 600 fits. The paper's other six algorithms (SNPE, SNLE, NRE, SNRE, REJ-ABC, SMC-ABC) are out of scope: `neuralsbi` has no ratio estimator and no ABC, and its sequential support is not the multi-round SNPE-C/SNLE-A the paper measured.

Nothing here is re-derived. Observations, true parameters and reference posteriors come straight out of the `sbibm` checkout, and the target numbers come straight out of `main_paper.csv` in `sbi-benchmark/results`. `neuralsbi` is scored against exactly the material the paper scored against.

## What you need

- R with `neuralsbi` installed (`R CMD INSTALL .` from the package root).
- `torch` with libtorch: `install.packages("torch"); torch::install_torch()`. The estimators the paper used are NSF for NPE and MAF for NLE, both torch-backed. There is no torch-free path through a real run.
- `deSolve`, for the two ODE tasks: `install.packages("deSolve")`. The other eight tasks do not need it.
- `git`, to fetch the two upstream repositories. No Python is needed at any point, including for the pickled tensors that three of the tasks depend on.
- Disk: about 40 MB for the upstream checkouts, plus a few hundred MB if you keep posterior samples for the full grid (`--save-samples FALSE` if you would rather not).

## Running it

```sh
cd dev/benchmarks

# 1. Fetch sbibm and the published results, and check everything is readable.
Rscript 00_setup.R

# 2. Check the harness before spending hours on it. Minutes.
Rscript smoke_test.R

# 3. Run something small first.
Rscript 01_run_benchmark.R --tasks two_moons --algorithms npe \
    --budgets 1000 --observations 1

# 4. Run the grid. Resumable: rerun after an interruption and it continues.
Rscript 02_run_all.R

# 5. Read the verdict.
Rscript 03_report.R
```

Every argument takes a comma-separated list, so you can carve the grid however you like:

```sh
Rscript 01_run_benchmark.R --tasks gaussian_linear,two_moons,slcp \
    --algorithms npe,nle --budgets 1000,10000 --observations 1,2,3
```

`02_run_all.R` runs one task per subprocess, cheapest first, with `--skip-existing` set. A task that fails is logged to `results/logs/<task>.log` and the run continues. Use `--dry-run` to print the commands it would issue.

If you already have checkouts of `sbibm` or `sbi-benchmark/results`, point `SBIBM_PATH` and `SBIBM_RESULTS_PATH` at them and skip `00_setup.R`. Set `NEURALSBI_BENCH_OUT` to write results somewhere other than `results/`.

## Reading the report

`03_report.R` prints one line per (task, algorithm, budget):

```
                   task alg  sims obs  ours     range paper  delta verdict
        gaussian_linear NPE  10^3  10 0.691 0.66-0.72 0.686 +0.005   MATCH
              two_moons NPE  10^3  10 0.702 0.66-0.77 0.715 -0.013   MATCH
```

The metric is C2ST: the cross-validated accuracy of a classifier trained to tell the approximate posterior from the reference posterior, on 10k samples from each. 0.5 means indistinguishable, 1.0 means perfectly separable, so lower is better. `ours` and `paper` are medians over the observations you ran, and the paper's median is taken over *those same observations*, so a partial run is still a fair comparison.

The verdict is `MATCH` when we are within the tolerance (0.05 by default, `--tolerance` to change it), `BETTER` when we are more than a tolerance below the paper, `WORSE` when we are more than a tolerance above it. Being better than the paper is a pass. The script ends with a count of reproduced cells and, if any cell is worse, prints those cells on their own. `--detail` adds the per-observation table.

Two files are written: `results/comparison.csv` (the per-cell table) and `results/report.md` (the same thing in markdown).

The `range` column is the min and max over the observations you ran. Read it before you read the verdict: on two_moons at 10³ simulations the per-observation C2ST spans roughly 0.65 to 0.77 in both implementations, so a median from three observations carries several hundredths of noise on its own.

Why 0.05: C2ST on 10k versus 10k samples has a run-to-run standard deviation of a few thousandths, but the training seed, the simulation seed and the observation all move the fitted posterior, and the paper reports a single run per cell. Differences under 0.05 do not tell you anything about the two implementations. Differences well over it do, and those are the cells worth opening up. If a cell lands near the boundary, rerun it with a different `--seed` before concluding anything; the results file keys on the cell, not the seed, so use a separate `NEURALSBI_BENCH_OUT` for each seed.

## How faithful the reimplementation is

The hyperparameters come from `config/algorithm/{npe,nle}.yaml` in the results repository, which is what produced `main_paper.csv`:

| | paper | here |
| --- | --- | --- |
| NPE estimator | `nsf`, 50 hidden features | `density_estimator = "nsf"`, `hidden = c(50, 50)` |
| NLE estimator | `maf`, 50 hidden features | `density_estimator = "maf"`, `hidden = c(50, 50)` |
| rounds | 1 | single round, which is what `npe()`/`nle()` do |
| z-scoring | `theta` and `x` | `standardize = TRUE` (the default) |
| training batch | 10000 | `batch_size = min(10000, n)` |
| NLE sampler | vectorized slice, 100 chains, thin 10, 100 warmup | `posterior(sampler = "slice", n_chains = 100, thin = 10, warmup = 100)` |
| posterior samples | 10000 | 10000 |

The rest of the training loop (5 transforms, 10 spline bins, tail bound 3, Adam at 5e-4, 10% validation split, early stopping after 20 epochs) is shared: `sbi`'s defaults and `neuralsbi`'s defaults agree there, which is not an accident.

The tasks in `R/tasks.R` follow `sbibm/tasks/<name>/task.py` term by term, including the frozen constants: the Bernoulli GLM's input stimulus and design matrix, and the SLCP distractor noise distribution, are read out of sbibm's own PyTorch files by `R/pt_io.R` rather than regenerated, because they came out of numpy's RNG and cannot be reproduced in R. The smoke test verifies the GLM path end to end by recomputing sbibm's published summary statistics from its published raw spike trains, and it verifies every simulator by checking that each shipped observation looks like a draw from our simulator at the parameters that generated it.

Four things genuinely differ, and they are worth knowing when a cell comes out `WORSE`:

1. **SIR and Lotka-Volterra integrate their ODEs with a different solver.** sbibm calls Julia's DifferentialEquations one parameter set at a time; we call `deSolve::ode()` the same way, which is `lsoda`: adaptive step, automatic stiff/non-stiff switching, `rtol` and `atol` both 1e-6 by default against Julia's 1e-3/1e-6. Output times are sbibm's `saveat` grid, subsampled exactly as sbibm subsamples it, and a solve that fails or returns anything non-finite becomes a NaN row, which is also what sbibm does. This is the slowest part of those two tasks at about 2 ms per solve, so set `NEURALSBI_BENCH_CORES` to fan the solves across cores; only the deterministic solves are parallelised, so the draws for a given seed do not depend on the core count.
2. **NLE runs in the original parameter space.** sbibm transforms parameters to an unconstrained space before training the likelihood and sampling it (`automatic_transforms_enabled: true`). `neuralsbi` keeps the parameters as they are and lets the prior's support cut the potential. This matters most for the log-normal priors of SIR and Lotka-Volterra.
3. **The two flow implementations are not the same code.** `neuralsbi`'s NSF is autoregressive; `sbi`'s is coupling-based. The MAFs are closer but still independent implementations.
4. **The C2ST classifier is a reimplementation.** `R/c2st.R` mirrors sbibm's scikit-learn `MLPClassifier` recipe (two hidden layers of 10 x dim, ReLU, Adam, 5-fold cross-validation, accuracy, z-scored by the reference sample) down to the initialisation scheme and the convergence rule, but sklearn's RNG cannot be reproduced outside Python. The epoch cap is 1000 rather than sbibm's 10000, for speed; raise it with `--c2st-max-epochs`. That cap biases in the safe direction, since a classifier cut short reports a *higher* accuracy on posterior samples this close together, which can only make `neuralsbi` look further from the reference than it is. Note also that `neuralsbi::c2st()` is a *different* test, a logistic regression, and its numbers are not comparable to the paper's; the scripts do not use it.

## What it costs

A single NPE cell at 10³ simulations is seconds to a couple of minutes. At 10⁵ the simulation itself dominates for the ODE tasks (a few minutes each) and training dominates for the rest. NLE adds an MCMC run per cell, and 100 chains of a slice sampler in R is not fast. Scoring is not free either, and on the wider tasks it dominates: the C2ST classifier on 10k against 10k samples costs roughly half a second per epoch per fold in ten dimensions, and it runs to the epoch cap rather than converging early, so budget on the order of twenty minutes per cell for the ten-parameter tasks. The two-parameter tasks are a minute or less. Lower `--c2st-max-epochs` if that is the binding constraint, and read the note on it below before you do. The full grid is a multi-day run on one machine. Realistic ways to cut it down:

- `--observations 1,2,3` instead of all ten. Three observations still give you a median with useful resolution.
- `--budgets 1000,10000`, leaving 10⁵ for a final confirmation.
- `--algorithms npe` first. NPE is the cheaper half by a wide margin.

Cells are independent, so several `01_run_benchmark.R` processes on disjoint task lists work fine. They all append to the same `results/metrics.csv`, so do not run two processes over the *same* cells at once.

## Files

```
00_setup.R            fetch sbibm + sbi-benchmark/results, verify readability
01_run_benchmark.R    run and score cells; appends to results/metrics.csv
02_run_all.R          drive the grid, one subprocess per task, resumable
03_report.R           compare against the paper, write comparison.csv + report.md
smoke_test.R          fast structural checks; run this before anything long

R/utils.R             argument parsing, paths, logging
R/pt_io.R             read tensors out of legacy PyTorch files, without Python
R/sbibm_data.R        observations, reference posteriors, published results
R/tasks.R             the ten sbibm tasks, vectorized
R/c2st.R              sbibm's C2ST, reimplemented to match
R/runner.R            one cell end to end, plus result bookkeeping

results/metrics.csv   one row per cell: C2ST, timings, acceptance rate
results/samples/      posterior draws, so a cell can be re-scored without refitting
results/comparison.csv, results/report.md   the verdict
```

`external/` and `results/` are git-ignored.

## Provenance

- Benchmark and reference material: <https://github.com/sbi-benchmark/sbibm>
- Published results: <https://github.com/sbi-benchmark/results>, `benchmarking_sbi/results/main_paper.csv`
- Paper: Lueckmann, J.-M., Boelts, J., Greenberg, D., Gonçalves, P. and Macke, J. (2021). Benchmarking Simulation-Based Inference. *AISTATS*, PMLR 130:343-351.
