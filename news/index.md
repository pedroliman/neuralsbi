# Changelog

## neuralsbi 0.4.4

- Fixed: an NPE posterior given a multi-row observation now warns.
  [`resolve_x()`](https://neuralsbi.pedrodelima.com/reference/resolve_x.md)
  kept row 1 and said nothing about the rest, which is right for NPE and
  invisible to the user. The trap is the asymmetry with
  [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md), where
  the rows of `x_obs` are independent observations that the
  log-likelihood sums over: the same matrix means “200 observations” to
  [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) and “the
  first observation” to
  [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md), so
  moving a working call from one to the other silently dropped all but
  one data point. The warning names the row count and points at
  [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md). It
  stays a warning because taking row 1 of a simulation matrix is a
  reasonable thing to ask for, and
  [`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md)/[`tarp()`](https://neuralsbi.pedrodelima.com/reference/tarp.md)
  pass single rows anyway.

## neuralsbi 0.4.3

- **New: Neural Likelihood Estimation, via
  [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md).** Where
  [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) learns
  the posterior directly,
  [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) learns a
  surrogate likelihood . The reason to want that is repeated
  observations. An NPE fit is trained for one fixed data dimension, so
  conditioning on independent trials means retraining for every or
  compressing to summary statistics; NLE learns the density of a single
  trial, so the log-likelihood of of them is a sum and is free at
  inference time. The cost is that the posterior is no longer a forward
  pass. For one fixed observation with high-dimensional data,
  [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) is still
  the better choice, and
  [`?nle`](https://neuralsbi.pedrodelima.com/reference/nle.md) says so.
- `log_lik(fit, theta, x)` evaluates the surrogate likelihood, summing
  over the rows of `x` as independent observations, and
  `likelihood_fn(fit, x_obs)` returns it as a plain vectorized
  `function(theta)`. That closure is the point of contact with the rest
  of R: it goes straight into
  [`optim()`](https://rdrr.io/r/stats/optim.html), an MCMC package, an
  importance sampler, or a profile likelihood, with nothing downstream
  needing to know about `neuralsbi`.
- [`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md)
  on an [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md)
  fit returns an MCMC-backed posterior and gains `sampler`, `n_chains`,
  `warmup`, `thin` and `init_strategy`. The default sampler is a
  vectorized univariate slice sampler, matching Python `sbi`’s default:
  nothing to tune, no dependency, and bounded priors need no special
  handling. Its vectorization runs across chains, so one step costs one
  batched forward pass rather than `n_chains` separate ones. The slice
  width adapts to the target during warmup, which matters because a
  width inherited from the prior is badly wrong once many observations
  have concentrated the posterior, and every unit of mismatch is paid
  for in wasted density evaluations. Draws carry split-Rhat and bulk
  ESS, and are cached on the posterior so
  [`summary()`](https://rdrr.io/r/base/summary.html) and repeat
  [`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md)
  calls do not re-run a chain. `n_chains` defaults to 20 under `"slice"`
  and 4 under `"stan"`, because a slice chain rides along in a batch
  someone else is already paying for and a Stan chain is a process with
  its own warmup.
- `thin` defaults to 2. With the adapted width, `thin = 2` already gives
  a bulk ESS of about 96% of the retained draws on a Gaussian target, so
  thinning harder buys very little for what it costs, and the reported
  ESS is there to tell you when to raise it. Python `sbi` thins by 1 (it
  thinned by 10 up to v0.21), so this sits between the two.
- The i.i.d. sum takes a shortcut where the estimator allows one. An MDN
  maps `theta` to a Gaussian mixture over `x` and never sees `x`, so for
  `n` observations the network runs once and all `n` densities come off
  the same mixture; the linear-Gaussian baseline behaves the same way. A
  flow’s transforms depend on `x` too, so it has to run `n` times. With
  a few thousand observations that is the difference between seconds and
  minutes per MCMC step, and it is worth weighing when choosing an
  estimator for repeated data.
- Sampling an NLE posterior is about five times faster than it was when
  the feature first worked, and samples exactly the same thing. Four
  changes get it there. The slice sampler now carries several of a
  chain’s future moves in one call: where an interval edge goes next is
  a fixed step, and which point a chain tries after a rejection depends
  only on that point’s side of the current value and a fresh uniform, so
  neither needs a density and both can be computed in advance. No call
  is ever made wider than the one that opened the coordinate, so this
  spends batch width that was already being paid for rather than adding
  any. Everything the observation alone decides – standardizing it, its
  Jacobian, turning it into a tensor – is now settled when the posterior
  is built instead of at every step. The summed log-likelihood is
  reduced where the densities are produced, so the `n_theta x n_obs`
  matrix, the largest object in the loop, is never built. And an MDN’s
  chunking is sized by the pair count rather than by pairs times
  components, which had been splitting a 5000-observation call ten ways
  for a 400 KB intermediate. On the four-parameter g-and-k model in
  [`vignette("neural-likelihood")`](https://neuralsbi.pedrodelima.com/articles/neural-likelihood.md),
  2000 draws from 20 chains went from 117s to 25s at 500 observations
  and from 342s to 68s at 5000, with split-Rhat and bulk ESS unchanged.
- An MDN likelihood is replayed as TorchScript once a run is clearly a
  loop. Every torch operation crosses from R into libtorch, and at MCMC
  batch sizes that crossing costs more than the arithmetic behind it:
  roughly 0.2 ms each and thirty per evaluation, against a few hundred
  microseconds of real work.
  [`torch::jit_trace()`](https://torch.mlverse.org/docs/reference/jit_trace.html)
  records the same code and replays it in one crossing. It is a
  shortcut, not a path. Nothing is recorded until an evaluator has been
  called a few times, so a single
  [`log_lik()`](https://neuralsbi.pedrodelima.com/reference/log_lik.md)
  never pays for a compiler; a trace fixes the shapes it saw, so there
  is one per parameter-row count and each is checked against the eager
  result before anything uses it; and tracing or checking failing just
  leaves the eager path in place. `options(neuralsbi.jit = FALSE)` turns
  it off.
- [`log_prob()`](https://neuralsbi.pedrodelima.com/reference/log_prob.md)
  on an NLE posterior returns the **unnormalized** log posterior. The
  evidence is not available, so `normalize` is ignored with a warning
  rather than returning a number that looks normalized and is not.
- **New:
  [`stan_code()`](https://neuralsbi.pedrodelima.com/reference/stan_export.md),
  [`stan_data()`](https://neuralsbi.pedrodelima.com/reference/stan_export.md)
  and
  [`write_stan_model()`](https://neuralsbi.pedrodelima.com/reference/stan_export.md)
  export a fitted likelihood as Stan source.** The generated `functions`
  block recomputes in Stan’s own language with the trained weights
  passed as data, so Stan differentiates it and NUTS gets exact
  gradients, with nothing linked against `torch` at run time. The result
  is an ordinary Stan function of `theta`, which is what makes it worth
  having: the surrogate stops being the whole model and becomes one term
  in a model you write, next to a hierarchical prior, covariates, or a
  second data source whose likelihood you do know. `"mdn"`, `"maf"` and
  `"linear_gaussian"` are supported; `"nsf"` is refused with a message
  naming the alternatives. `posterior(fit, x_obs, sampler = "stan")`
  runs the generated model through `cmdstanr` or `rstan` and returns
  draws like any other path.
- [`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md),
  [`tarp()`](https://neuralsbi.pedrodelima.com/reference/tarp.md) and
  [`posterior_predictive()`](https://neuralsbi.pedrodelima.com/reference/posterior_predictive.md)
  accept an
  [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) fit, and
  [`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md)/[`tarp()`](https://neuralsbi.pedrodelima.com/reference/tarp.md)
  forward `...` to
  [`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md)
  so the MCMC controls reach it. Every SBC trial is a separate MCMC run,
  so start small.
- [`save_npe()`](https://neuralsbi.pedrodelima.com/reference/save_npe.md)/[`load_npe()`](https://neuralsbi.pedrodelima.com/reference/save_npe.md)
  handle [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md)
  fits too, with
  [`save_nle()`](https://neuralsbi.pedrodelima.com/reference/save_npe.md)/[`load_nle()`](https://neuralsbi.pedrodelima.com/reference/save_npe.md)
  as aliases.
- Internally,
  [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) and
  [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) now
  share
  [`prepare_simulations()`](https://neuralsbi.pedrodelima.com/reference/prepare_simulations.md)
  instead of each carrying its own copy of the
  simulate/coerce/drop/standardize preamble.
- Fixed:
  [`sample_posterior()`](https://neuralsbi.pedrodelima.com/reference/sample_posterior.md)
  now goes through the
  [`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md)
  generic instead of calling
  [`sample.nsbi_posterior()`](https://neuralsbi.pedrodelima.com/reference/sample.nsbi_posterior.md)
  directly. An
  [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md)
  posterior inherits `nsbi_posterior`, so the old call ran the NPE
  forward-pass sampler against an estimator whose target and condition
  are swapped. That errored with `non-conformable arguments` when
  `dim_theta` and `dim_x` differ, and returned draws from the wrong
  distribution without complaint when they happen to match.
- Fixed:
  [`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md)
  on an [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md)
  fit now validates `n_chains`, `warmup` and `thin` instead of coercing
  them with [`as.integer()`](https://rdrr.io/r/base/integer.html).
  `thin = 0` ran the warmup and kept nothing, so the slice sampler
  returned its zero-initialized array and every draw was 0, with no
  error and an Rhat computed on those zeros. Values below the bound,
  fractional values and `NA` are now errors.
- Fixed: the non-finite check now covers `theta` as well as `x`. A row
  was dropped only when its simulator output was non-finite, so an `NA`
  among pre-computed parameters passed straight through to the
  estimator: `npe(prior, theta = theta, x = x)` with one missing
  parameter failed inside [`chol()`](https://rdrr.io/r/base/chol.html)
  with “the leading minor of order 1 is not positive”, and on the torch
  path it became “Training failed: no restart produced a finite
  validation loss”, which blames training for a bad input. Such rows are
  now dropped like any other failure, and the warning says which side
  was non-finite (“parameters”, “output”, or “parameters or output”).
- Fixed:
  [`npe_sequential()`](https://neuralsbi.pedrodelima.com/reference/npe_sequential.md)
  now checks `x_obs` against the simulator’s output width, at the end of
  round 1, which is the first moment that width is known. An `x_obs` of
  the wrong length was reshaped into several rows and only the first was
  targeted, so the run truncated its proposals around a value the caller
  never gave while [`print()`](https://rdrr.io/r/base/print.html)
  reported all of them as targeted. A TSNPE fit is not amortized, so
  targeting the wrong observation is the whole failure. A missing,
  `NULL`, non-numeric, `NA` or multi-row `x_obs` is an error as well,
  and `n_rounds` must be a single integer of at least 1: `n_rounds = 0`
  skipped the loop and returned a bare list classed `nsbi_snpe`.
- Fixed: [`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md)
  and [`tarp()`](https://neuralsbi.pedrodelima.com/reference/tarp.md)
  now error when a trial’s posterior returns fewer draws than
  `n_posterior_samples` instead of scoring it anyway.
  [`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md)
  comes back short when a bounded prior and an estimator that leaks mass
  outside it defeat rejection sampling, and it only warns.
  [`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md) then
  binned those ranks against `n_posterior_samples` while they were drawn
  from a smaller set, which compressed every rank toward zero: the
  chi-square uniformity test rejected and
  [`expected_coverage()`](https://neuralsbi.pedrodelima.com/reference/expected_coverage.md)
  fell below the diagonal, so lost draws were reported as a
  miscalibrated posterior. Rescaling the short trial on its own is no
  better, since it is then scored on a different resolution from the
  rest, so both functions stop and name the trial and the shortfall.
- Fixed: [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md)
  and [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) now
  match `density_estimator` against the allowed names before the
  simulator runs. The check lived in `fit_density_estimator()`, which is
  reached only after the simulations are in hand, so
  `density_estimator = "mfa"` spent the whole budget and then failed on
  the typo. For the expensive simulators this package is for, that is
  minutes to hours. The same hoist fixes the embedding warning, which
  compared the unmatched value:
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html) accepts
  abbreviations, so `density_estimator = "linear"` selected
  `linear_gaussian` and dropped `embedding_net` without saying so.
  `fit_density_estimator()` keeps its own
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html) because it is
  also reachable directly.
- Fixed: a simulator whose formals match some parameter names but not
  all now warns. The choice between the two simulator signatures was all
  or nothing and the fallback was silent, so
  `prior_uniform(c(a = -3, b = -1), c(a = 3, b = 1))` with
  `function(mu, ls = 0)` sent the whole two-parameter vector to `mu`,
  left `ls` at its default, and trained on the result without complaint.
  One typo in a prior name was enough. The warning names the parameters
  that found no formal and points at
  [`?nsbi_simulator`](https://neuralsbi.pedrodelima.com/reference/nsbi_simulator.md).
  It is a warning and not an error because a vector-signature simulator
  whose first argument happens to carry a parameter’s name is legal.

## neuralsbi 0.4.2

- [`pairplot()`](https://neuralsbi.pedrodelima.com/reference/pairplot.md)’s
  lower triangle now shows highest-density regions (via the new
  `ggdensity` Suggests, `geom_hdr()`, at the 50/80/95/99% probability
  levels) instead of a raw scatter of draws. A scatter of 10,000+ points
  overplots into an undifferentiated blob at the resolution most
  posteriors are viewed at; nested HDR contours show where the mass
  actually concentrates. The `truth` cross-hair markers are unchanged.
  `col` now sets the region fill colour; `alpha` applies only to the
  diagonal marginal densities, since the lower triangle shades itself by
  probability level.

## neuralsbi 0.4.1

- **Breaking: the simulator is now called once per parameter set and
  returns one simulated observation.** Most models a researcher already
  has work that way – an ODE solve, an agent-based model, a call out to
  `pomp` or `deSolve` – so meeting the old contract meant wrapping the
  model in an apply loop and thinking about column-major indexing before
  anything ran. Two signatures are accepted, decided from
  `formals(simulator)`: when every parameter name appears among the
  formals the parameters arrive by name, one scalar each
  (`function(mu, sigma) ...`); otherwise the whole named parameter
  vector goes to the first argument (`function(theta) ...`). A simulator
  returns a numeric vector of length `d`, a scalar, or a one-row matrix
  or data frame, and names on that output become the outcome names. See
  [`?nsbi_simulator`](https://neuralsbi.pedrodelima.com/reference/nsbi_simulator.md)
  for the contract and the migration examples.
- This also fixes a correctness bug introduced in 0.4.0. Chunking meant
  a simulator written for a single parameter set could return a vector
  whose length happened to match the chunk’s row count, and that vector
  was accepted as one column of outcomes for several draws.
  `simulate_for_sbi(function(theta) c(theta[1] * 2, theta[2] * 2), prior, 100)`
  returned a 100 x 1 matrix of nonsense with no error and no warning; it
  now returns the 100 x 2 matrix it should.
- [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md),
  [`simulate_for_sbi()`](https://neuralsbi.pedrodelima.com/reference/simulate_for_sbi.md),
  [`npe_sequential()`](https://neuralsbi.pedrodelima.com/reference/npe_sequential.md),
  [`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md),
  [`tarp()`](https://neuralsbi.pedrodelima.com/reference/tarp.md) and
  [`posterior_predictive()`](https://neuralsbi.pedrodelima.com/reference/posterior_predictive.md)
  gain `sim_args`, a named list forwarded to every simulator call.
  Observed data, a time grid, a fixed population size, a design matrix
  or a solver tolerance travel from the call site instead of being
  captured in a closure – which also keeps them out of what gets
  serialized to a `future` worker. A list rather than `...` because
  every one of
  [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md)’s
  formals sits before `...`, so R’s partial matching would silently
  capture `x`, `theta`, `n` or `seed`.
- Simulations whose output contains `NA`, `NaN` or an infinite value are
  dropped, together with their parameters, with one warning per run
  reporting the count and the rate. A single non-finite value would
  otherwise poison the training loss and surface much later as a `NaN`
  validation loss. The count is recorded on the fit and shown by
  [`print()`](https://rdrr.io/r/base/print.html). Nothing left is an
  error. In
  [`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md) and
  [`tarp()`](https://neuralsbi.pedrodelima.com/reference/tarp.md) a
  failed draw removes the whole trial; in
  [`posterior_predictive()`](https://neuralsbi.pedrodelima.com/reference/posterior_predictive.md)
  it reduces the number of predictive draws. Pre-computed `theta`/`x`
  passed to
  [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) are
  checked the same way. Note that dropping conditions on the simulator
  having succeeded: when failure depends on the parameters, the fit
  targets the posterior *given success*, so a high drop rate is a
  modelling signal.
- Random-number streams are now per simulation rather than per chunk, so
  a given seed produces the same simulations whatever the `future` plan
  and whatever the worker count. The 0.4.0 guarantee held only at a
  fixed chunk size.
- **`chunk_size` is gone**, from
  [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md),
  [`simulate_for_sbi()`](https://neuralsbi.pedrodelima.com/reference/simulate_for_sbi.md),
  [`npe_sequential()`](https://neuralsbi.pedrodelima.com/reference/npe_sequential.md),
  [`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md),
  [`tarp()`](https://neuralsbi.pedrodelima.com/reference/tarp.md) and
  [`posterior_predictive()`](https://neuralsbi.pedrodelima.com/reference/posterior_predictive.md),
  along with `options(neuralsbi.chunks)`. Chunking existed in 0.4.0 to
  make results reproducible across backends: the split had to depend on
  `n` alone, which meant users had to know about it and could change
  their draws by changing it. Per-simulation RNG streams give that
  guarantee outright, so what is left is ordinary parallel scheduling.
  Batches are now sized from the worker count, cannot affect a result,
  and are not a setting. Running sequentially there are no batches at
  all, just a loop.
- NSF’s `n_bins` and `tail_bound` are explicit
  [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md)
  arguments;
  [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) no
  longer takes `...`.
- New
  [`save_npe()`](https://neuralsbi.pedrodelima.com/reference/save_npe.md)
  /
  [`load_npe()`](https://neuralsbi.pedrodelima.com/reference/save_npe.md).
  A fit whose estimator is `"maf"`, `"mdn"` or `"nsf"` holds a torch
  module, which is an external pointer:
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) writes the pointer,
  the file reloads without complaint, and the first call that touches
  the network fails with `external pointer is not valid`.
  [`save_npe()`](https://neuralsbi.pedrodelima.com/reference/save_npe.md)
  writes the weights with
  [`torch::torch_save()`](https://torch.mlverse.org/docs/reference/torch_save.html)
  and everything else as ordinary R objects, into one `.rds`;
  [`load_npe()`](https://neuralsbi.pedrodelima.com/reference/save_npe.md)
  rebuilds the network from the recorded architecture and restores them.
  An overnight fit can now be reloaded the next morning, which is what
  amortization was for.
  [`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md)
  and [`print()`](https://rdrr.io/r/base/print.html) also detect a fit
  that came back from [`readRDS()`](https://rdrr.io/r/base/readRDS.html)
  and say so, instead of failing later with a torch error.

## neuralsbi 0.4.0

- The simulator can now run in parallel. Declare a `future` plan –
  [`library(future); plan(multisession)`](https://future.futureverse.org)
  – and every function that calls a simulator
  ([`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md),
  [`simulate_for_sbi()`](https://neuralsbi.pedrodelima.com/reference/simulate_for_sbi.md),
  [`npe_sequential()`](https://neuralsbi.pedrodelima.com/reference/npe_sequential.md),
  [`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md),
  [`tarp()`](https://neuralsbi.pedrodelima.com/reference/tarp.md),
  [`posterior_predictive()`](https://neuralsbi.pedrodelima.com/reference/posterior_predictive.md))
  spreads the work across workers. There is no new argument to pass and
  no parallel variant to call: with no plan declared everything runs
  sequentially as before, and `neuralsbi` mentions the two lines above
  once per session (`options(neuralsbi.parallel_hint = FALSE)` to
  silence it). Each chunk of parameters draws from its own L’Ecuyer-CMRG
  stream, so a given [`set.seed()`](https://rdrr.io/r/base/Random.html)
  produces the same simulations sequentially and on any number of
  workers. See
  [`?nsbi_parallel`](https://neuralsbi.pedrodelima.com/reference/nsbi_parallel.md).

- Long-running work now reports progress with an ETA – simulation and
  neural training alike, one progress step per simulation and per
  training epoch. With `progressr` installed, `neuralsbi` emits standard
  progressr updates, so
  [`progressr::handlers()`](https://progressr.futureverse.org/reference/handlers.html)
  and
  [`with_progress()`](https://progressr.futureverse.org/reference/with_progress.html)
  control reporting; without it, a built-in bar needing no extra
  packages does the job. The training bar targets the epoch at which
  early stopping would fire and revises that target as the validation
  loss improves. See
  [`?nsbi_progress`](https://neuralsbi.pedrodelima.com/reference/nsbi_progress.md).

- [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md),
  [`simulate_for_sbi()`](https://neuralsbi.pedrodelima.com/reference/simulate_for_sbi.md),
  [`npe_sequential()`](https://neuralsbi.pedrodelima.com/reference/npe_sequential.md),
  [`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md),
  [`tarp()`](https://neuralsbi.pedrodelima.com/reference/tarp.md), and
  [`posterior_predictive()`](https://neuralsbi.pedrodelima.com/reference/posterior_predictive.md)
  gain a `chunk_size` argument controlling how many parameter rows go to
  the simulator per call. The default splits a run into about 64 chunks;
  because the split depends only on the number of simulations, results
  do not change with the number of workers. A simulator whose output
  must be produced in one call can set `chunk_size` to the full
  simulation budget.

- `future` and `progressr` are `Suggests`, not dependencies; `parallel`
  (base R) moves into `Imports`. \# neuralsbi 0.3.7

- The test suite now skips its plotting tests when `ggplot2`/`GGally`
  are not installed, instead of failing. Both are `Suggests`, so
  `R CMD check` under `_R_CHECK_FORCE_SUGGESTS_=false` – the
  configuration CRAN uses on a machine without them – previously hit 5
  errors from
  [`require_ggplot2()`](https://neuralsbi.pedrodelima.com/reference/require_ggplot2.md).
  The new `skip_if_no_ggplot2()`/`skip_if_no_ggally()` helpers mirror
  the `skip_if_no_torch()` contract already used for the neural tests,
  so the suite runs everywhere.

- [`vignette("sir-time-varying-beta")`](https://neuralsbi.pedrodelima.com/articles/sir-time-varying-beta.md)
  reworks the epidemic model so that its posterior predictive tracks the
  observed case peaks instead of overshooting them. The introduction day
  and seed size are now inferred per state rather than fixed at two
  infections on 2020-01-21 (which left the epidemic’s phase to
  demographic noise: replicate simulations at one fixed parameter
  differed by a factor of several thousand at the peak bin), the
  ascertainment prior is pinned by spring-2020 seroprevalence (case
  counts alone identify only the product of ascertainment and incidence,
  and the unconstrained fit slides toward vast, barely-ascertained
  epidemics), and the regime-duration parameterization drops a
  coordinate that the simulator’s rescaling had left unidentified. The
  article reports the effective reproduction number *R*_(e)(t) = beta(t)
  S(t) / (N gamma) computed from the simulator’s own susceptible
  trajectories, adds a prior-predictive check before training, and
  summarizes the posterior predictive by its median rather than its
  mean.

## neuralsbi 0.3.6

- Parameters and outcomes can now be named. Name
  [`prior_uniform()`](https://neuralsbi.pedrodelima.com/reference/prior_uniform.md)’s
  `low`/`high` or
  [`prior_normal()`](https://neuralsbi.pedrodelima.com/reference/prior_normal.md)’s
  `mean` (e.g. `c(beta = 0, gamma = 0)`), or attach
  [`colnames()`](https://rdrr.io/r/base/colnames.html) to a simulator’s
  output, and those names carry through
  [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md),
  [`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md),
  [`map_estimate()`](https://neuralsbi.pedrodelima.com/reference/map_estimate.md),
  [`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md),
  [`expected_coverage()`](https://neuralsbi.pedrodelima.com/reference/expected_coverage.md),
  and
  [`posterior_predictive()`](https://neuralsbi.pedrodelima.com/reference/posterior_predictive.md)
  without any extra arguments.
  [`plot_sbc()`](https://neuralsbi.pedrodelima.com/reference/plot_sbc.md),
  [`plot_coverage()`](https://neuralsbi.pedrodelima.com/reference/plot_coverage.md),
  [`pairplot()`](https://neuralsbi.pedrodelima.com/reference/pairplot.md),
  and
  [`plot_posterior_predictive()`](https://neuralsbi.pedrodelima.com/reference/plot_posterior_predictive.md)
  use the names for titles, legends, and facet labels; a name that
  happens to be valid R syntax (`"beta[1]"`, `"rho"`, `"sigma^2"`)
  renders as its plotmath symbol (Greek letters, sub/superscripts)
  instead of literal text.

## neuralsbi 0.3.5

- The SIR case study becomes a head-to-head comparison with the `pomp`
  package,
  [`vignette("sir-epidemic")`](https://neuralsbi.pedrodelima.com/articles/sir-epidemic.md).
  Both methods fit the same stochastic SIR epidemic: `pomp` via
  particle-filter MCMC (`pmcmc`), `neuralsbi` via neural posterior
  estimation. The vignette contrasts what each needs from the model —
  `pomp` a measurement density, `neuralsbi` only a simulator — overlays
  the two posteriors, scores their agreement with a C2ST, and confirms
  the neural fit with SBC. The comparison is precomputed, so `pomp` is
  needed only to regenerate the article, not to build or check the
  package.

## neuralsbi 0.3.4

- Plotting is now built on `ggplot2` and
  [`GGally::ggpairs()`](https://ggobi.github.io/ggally/reference/ggpairs.html)
  instead of base graphics.
  [`pairplot()`](https://neuralsbi.pedrodelima.com/reference/pairplot.md),
  [`plot_sbc()`](https://neuralsbi.pedrodelima.com/reference/plot_sbc.md),
  [`plot_coverage()`](https://neuralsbi.pedrodelima.com/reference/plot_coverage.md),
  [`plot_tarp()`](https://neuralsbi.pedrodelima.com/reference/plot_tarp.md),
  and
  [`plot_posterior_predictive()`](https://neuralsbi.pedrodelima.com/reference/plot_posterior_predictive.md)
  keep their signatures
  ([`pairplot()`](https://neuralsbi.pedrodelima.com/reference/pairplot.md)
  gains an `alpha` argument) but now build and print a
  `ggplot`/`ggmatrix` object, returned invisibly for further
  customization. `ggplot2` and `GGally` move to `Suggests`, following
  the same graceful-degradation pattern as `torch`: call any plotting
  function without them installed and you get an informative error, not
  a crash.
- New vignette,
  [`vignette("intro-to-sbi")`](https://neuralsbi.pedrodelima.com/articles/intro-to-sbi.md):
  a short beginner tutorial covering the three ingredients (prior,
  simulator, observation), amortized training, and a first calibration
  check, using a g-and-k distribution simulator whose likelihood has no
  closed form.

## neuralsbi 0.3.2

- README is now generated from `README.Rmd`, so the usage example runs
  at render time and its output cannot drift from the code. The example
  is a plain linear regression: the posterior recovers the ground-truth
  coefficients, which the reader can check against ordinary least
  squares.

## neuralsbi 0.3.1

- CRAN resubmission fixes. Routed the
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods
  into the single `summaries` help topic (they had drifted into separate
  `.Rd` files with duplicated `\alias` entries, which also produced
  duplicate HTML anchors). Wrapped `theta_{<d}` in `\eqn{}` in the
  `made_masks` docs so the Rd no longer drops braces. Dropped the bare
  “NPE” acronym from the `DESCRIPTION` to avoid the spurious misspelling
  note.

## neuralsbi 0.3.0

- Defaults now match Python `sbi`, so a workflow reads the same in both
  packages and results can be cross-checked. Changes to
  [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md)
  defaults: the density estimator is now `"maf"` (was `"mdn"`); MDN
  mixture components default to 10 (was 5); NSF spline bins default to
  10 (was 8); the training batch size is 200 (was 100). `max_epochs` is
  raised to 2000 as a guard cap that early stopping (`patience = 20`)
  normally reaches first, mirroring `sbi`’s effectively-unbounded epoch
  budget. `lr`, `validation_fraction`, `patience`, `clip_grad_norm`,
  `n_transforms`, and `hidden` already matched. Pass any of these
  explicitly to recover the previous behavior.

- First CRAN submission. Dropped the development `.9000` version suffix,
  removed the redundant `Author`/`Maintainer` fields (now derived from
  `Authors@R`), and tidied the package title.

- Embedding networks (roadmap v0.4).
  [`embedding_mlp()`](https://neuralsbi.pedrodelima.com/reference/embedding_mlp.md)
  builds a learned summary network that maps raw observations to a
  low-dimensional feature vector; pass it to
  `npe(..., embedding_net = )` and the MDN, MAF, and NSF estimators
  condition on the features instead of the raw data, training the
  embedding jointly. The estimators still take raw `x` at the `de_*`
  boundary (`dim_x` is unchanged), so sampling and `log_prob` route
  through the embedding automatically. Ignored, with a warning, by
  `linear_gaussian`.

## neuralsbi 0.2.4.9000 (development)

- Vignettes now show real output. They are *precomputed*: each
  vignette’s evaluated source lives in `vignettes/<name>.Rmd.orig`, and
  `vignettes/precompute.R` bakes it into a static `vignettes/<name>.Rmd`
  (results, printed values, and figures inlined). CI and pkgdown
  re-render that static Markdown with no torch at build time, so the
  expensive neural training runs once, locally, instead of on every
  build. Re-run `Rscript vignettes/precompute.R` after editing any
  `.Rmd.orig`.
- Two-moons calibration study
  (`inst/benchmarks/two_moons_calibration.R`): SBC, expected coverage,
  and TARP for a two-moons NSF fit, with figures written to
  `docs/figures/` (roadmap milestone M2).

## neuralsbi 0.2.3.9000

- Package website built with pkgdown, deployed from CI to
  <https://pedroliman.github.io/neuralsbi/>.
- Four vignettes that build on each other: getting started, choosing a
  density estimator, checking the posterior, and the SIR case study
  (which now also demonstrates
  [`npe_sequential()`](https://neuralsbi.pedrodelima.com/reference/npe_sequential.md)).
  Removed a truncated duplicate of the SIR vignette.
- README rewritten to the standard terse form; authorship recorded in
  `DESCRIPTION` (Pedro Nascimento de Lima, with ORCID).

## neuralsbi 0.2.2.9000

- New
  [`npe_sequential()`](https://neuralsbi.pedrodelima.com/reference/npe_sequential.md):
  multi-round NPE targeting a single observation via truncated-prior
  proposals (TSNPE, Deistler et al. 2022). Each round truncates the
  prior to the highest-probability region of the current posterior and
  retrains on all accumulated simulations; the standard NPE loss stays
  valid, so no importance correction is needed. Returns an `nsbi_snpe`
  fit that works with
  [`posterior()`](https://neuralsbi.pedrodelima.com/reference/posterior.md),
  [`sample()`](https://neuralsbi.pedrodelima.com/reference/sample.md),
  and the diagnostics, but is only valid at the targeted `x_obs`.
  Verified against the analytic linear-Gaussian posterior.

## neuralsbi 0.2.1.9000

- New [`tarp()`](https://neuralsbi.pedrodelima.com/reference/tarp.md)
  diagnostic and
  [`plot_tarp()`](https://neuralsbi.pedrodelima.com/reference/plot_tarp.md)
  (Lemos et al. 2023): a *joint* expected-coverage test using random
  reference points, complementing the per-parameter
  [`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md) ranks.
  Detects posteriors with calibrated marginals but wrong correlation
  structure.
- New
  [`plot_posterior_predictive()`](https://neuralsbi.pedrodelima.com/reference/plot_posterior_predictive.md):
  marginal predictive histograms with the observation marked; returns
  the observation’s predictive quantiles.
- Leakage correction is now under test: with a bounded prior, the
  renormalized
  [`log_prob()`](https://neuralsbi.pedrodelima.com/reference/log_prob.md)
  integrates to one over the support and returns `-Inf` outside it
  (`test-posterior-normalization.R`).
- Fixed CI. `R CMD check` failed on three counts: the
  [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) example
  required libtorch (it now uses the torch-free `linear_gaussian`
  estimator and runs unconditionally), the hand-maintained
  `npe.Rd`/`fit_mdn.Rd` usage sections had drifted behind the code
  (missing `n_restarts`, `clip_grad_norm`, `n_transforms`, and the
  `"maf"`/`"nsf"` options), and `CLAUDE.md` was not in `.Rbuildignore`.
  The `test-torch` job also failed because torch 0.17 refuses a
  `TORCH_HOME` that does not exist; the workflow now creates it first.

## neuralsbi 0.2.0.9000

- Shared training engine for all neural estimators
  ([`train_conditional_de()`](https://neuralsbi.pedrodelima.com/reference/train_conditional_de.md)):
  best-of-n restarts, learning-rate decay on plateau, gradient clipping,
  per-epoch loss history.
- Masked Autoregressive Flow (`density_estimator = "maf"`) and Neural
  Spline Flow (`"nsf"`, autoregressive rational-quadratic splines) join
  the MDN and the closed-form `linear_gaussian` baseline.
- Benchmark tasks
  ([`task_gaussian_linear()`](https://neuralsbi.pedrodelima.com/reference/tasks.md),
  [`task_two_moons()`](https://neuralsbi.pedrodelima.com/reference/tasks.md),
  [`task_slcp()`](https://neuralsbi.pedrodelima.com/reference/tasks.md),
  [`task_sir()`](https://neuralsbi.pedrodelima.com/reference/tasks.md))
  shared between tests and the `inst/benchmarks/` head-to-head benchmark
  harness.
- [`summary()`](https://rdrr.io/r/base/summary.html) methods,
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) tidy
  accessor,
  [`plot_coverage()`](https://neuralsbi.pedrodelima.com/reference/plot_coverage.md).
- SIR applied case-study vignette.
- CI: `R CMD check` plus a `test-torch` job with cached libtorch.

## neuralsbi 0.1.0

- First pilot release: priors, single-round amortized
  [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md),
  `linear_gaussian` and MDN estimators, posterior sampling with leakage
  correction, SBC, expected coverage, C2ST, posterior-predictive checks,
  [`pairplot()`](https://neuralsbi.pedrodelima.com/reference/pairplot.md),
  [`plot_sbc()`](https://neuralsbi.pedrodelima.com/reference/plot_sbc.md).
