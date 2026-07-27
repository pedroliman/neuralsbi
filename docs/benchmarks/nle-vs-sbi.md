# NLE head-to-head against Python `sbi`

Run July 2026 with `neuralsbi` 0.4.3 and `sbi` 0.26.1 (torch 2.13, CPU).
Protocol and scripts: `inst/benchmarks/`, files `05`–`08`. Raw numbers:
`comparison_nle_maf.csv` beside this file.

## What was run

`task_gaussian_linear(dim = 5)`: prior `N(0, 0.1 I)`, likelihood
`x | theta ~ N(theta, 0.1 I)`. 10,000 simulations, one observation per
simulation, written once and read by both implementations, so neither sees a
draw the other did not.

Both trained a MAF likelihood on their own defaults (5 transforms, 50 hidden,
batch 200, lr 5e-4, 10% validation, early stopping after 20 epochs) and sampled
the resulting posterior with their own default MCMC (a vectorized slice sampler
in both), 5000 draws, at three observation-set sizes drawn from one parameter:
`n = 1`, `10`, `100`. Neither side was tuned. Tuning either one answers a
different question.

The point of the conjugate task is that the posterior stays closed-form for
every `n`: variance `1 / (1/v0 + n/s2)`, mean that variance times `sum(x) / s2`.
So both are scored against an exact answer, not against each other.

## Results

| `n_obs` | C2ST ours vs exact | C2ST sbi vs exact | C2ST ours vs sbi | max mean error, ours | max mean error, sbi | sd ratio, ours | sd ratio, sbi |
|---|---|---|---|---|---|---|---|
| 1   | 0.510 | 0.503 | 0.499 | 0.011 | 0.015 | 0.993 | 0.993 |
| 10  | 0.555 | 0.535 | 0.569 | 0.022 | 0.014 | 1.004 | 1.004 |
| 100 | 0.721 | 0.742 | 0.776 | 0.027 | 0.031 | 0.978 | 0.966 |

Read the C2ST next to the moments, not on its own: `c2st()` trains a logistic
regression, which sees a shift in location and is close to blind to a difference
in spread (`?c2st`). `sbi`'s trains an MLP.

**At one observation both implementations recover the analytic posterior and
are indistinguishable from each other** (C2ST 0.499, the coin-flip). Mean errors
of 0.011 and 0.015 sit against an exact posterior sd of 0.224, so both are
inside 7% of a standard deviation.

**At 100 observations both degrade, by the same amount, in the same way.** The
posterior widths stay right (sd ratios 0.98 and 0.97); what moves is the
location, by 0.027 and 0.031 against an exact posterior sd of 0.0315 -- about
0.85 and 0.99 standard deviations. The C2ST against each other (0.776) is no
better than either one's C2ST against the truth, which is what it looks like
when two estimates are wrong independently rather than one being wrong.

That degradation is the i.i.d. sum, not either implementation. A surrogate off
by `eps` nats per observation is off by `100 eps` over a hundred of them, while
the posterior sd contracts by a factor of ten. `vignette("neural-likelihood")`
runs into the same thing on the g-and-k model and spends a section on it; this
run is the controlled version, where the exact posterior is known at every `n`
and the effect can be measured rather than inferred.

## Reading the acceptance criterion

`inst/benchmarks/README.md` sets C2ST <= 0.60 against a reference posterior.
This run meets it at `n = 1` and `n = 10` and misses it at `n = 100`, for both
implementations. The criterion was written for NPE at a single observation, and
at `n = 100` it is no longer measuring an implementation: it is measuring the
surrogate's error per observation, multiplied by 100. What the criterion should
say for NLE, and what this run passes, is that the two implementations agree
with each other about as well as either agrees with the truth, and that both
recover the analytic posterior at `n = 1`.

## Cost

The two runs are not a clean timing comparison -- they shared four cores with
each other and with a vignette bake -- so the wall clocks (about 25 minutes for
`neuralsbi`, about 50 for `sbi`, covering training plus all three posteriors)
are indicative at best.

One number is not noise, because `sbi` reports it itself: **2 minutes 6 seconds
to generate 20 chain initializations**, before the first MCMC step. `sbi`'s
`resample` strategy draws `num_candidate_samples = 10_000` from the prior and
evaluates the potential on all of them, once per chain -- 200,000 evaluations at
20 chains. `neuralsbi` draws one pool of 1000, shares it across chains, and
takes a Gumbel top-k without replacement, which is 1000 evaluations and gives
distinct starting points rather than independent ones that may coincide. See
`mcmc_init()` in `R/mcmc.R`.

## Reproducing

```sh
cd inst/benchmarks
Rscript 05_generate_data_nle.R --dim 5 --n 10000 --seed 42
python  06_run_sbi_nle.py      --estimator maf --n_samples 5000 --seed 42
Rscript 07_run_neuralsbi_nle.R --estimator maf --n_samples 5000 --seed 42
Rscript 08_compare_nle.R       --estimator maf
```
