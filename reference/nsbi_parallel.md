# Running the simulator in parallel

Every `neuralsbi` function that calls a simulator –
[`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md),
[`simulate_for_sbi()`](https://neuralsbi.pedrodelima.com/reference/simulate_for_sbi.md),
[`npe_sequential()`](https://neuralsbi.pedrodelima.com/reference/npe_sequential.md),
[`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md),
[`tarp()`](https://neuralsbi.pedrodelima.com/reference/tarp.md),
[`posterior_predictive()`](https://neuralsbi.pedrodelima.com/reference/posterior_predictive.md)
– runs it through one execution path. That path is sequential by default
and parallel as soon as you declare a future plan:

## Details

    library(future)
    plan(multisession)          # or multicore, cluster, batchtools, ...
    fit <- npe(prior, simulator, n_simulations = 10000)

Nothing else changes: no argument to set, no variant of the function to
call. With no plan (or `plan(sequential)`) the simulations run in the
current process, and `neuralsbi` mentions the two lines above once per
session. Silence that hint with
`options(neuralsbi.parallel_hint = FALSE)`.

## Random numbers

Each *simulation* gets its own L'Ecuyer-CMRG random-number stream,
derived from the session's RNG state at the moment simulation starts. So
`set.seed(42)` before a fit (or `npe(..., seed = 42)`) gives the same
simulations whatever the plan and whatever the worker count: results
depend on the seed alone.

Because the simulations are drawn from separate streams, the simulator
no longer consumes the caller's RNG state directly; the state advances
by one draw per simulation phase.

## Scheduling

Under a multi-process plan the draws are dealt out to workers in
batches, because one `future` per simulation would spend longer shipping
the simulator to a worker than running it. Batch sizes follow the worker
count and are not a setting: they cannot change a result, since every
simulation has its own RNG stream. Running sequentially there are no
batches at all, just a loop.

Each batch ships the simulator, and everything its environment captures,
to a worker. A simulator closing over a large object therefore pays for
that object once per batch, which is a reason to pass it through
`sim_args` instead of capturing it.

## What is not parallelized

Training runs in the calling process. Torch tensors and modules are
external pointers that cannot be shipped to a worker, and libtorch
already uses multiple threads internally
([`torch::torch_set_num_threads()`](https://rdrr.io/pkg/torch/man/threads.html)).
The same goes for the posterior-sampling loops in
[`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md) and
[`tarp()`](https://neuralsbi.pedrodelima.com/reference/tarp.md) – those
call the trained network, not the simulator, and are reported with their
own progress bar.

## See also

[nsbi_progress](https://neuralsbi.pedrodelima.com/reference/nsbi_progress.md)
for controlling progress bars.
