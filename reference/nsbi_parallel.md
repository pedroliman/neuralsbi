# Running the simulator in parallel

Every `neuralsbi` function that calls a simulator –
[`npe()`](https://pedroliman.github.io/neuralsbi/reference/npe.md),
[`simulate_for_sbi()`](https://pedroliman.github.io/neuralsbi/reference/simulate_for_sbi.md),
[`npe_sequential()`](https://pedroliman.github.io/neuralsbi/reference/npe_sequential.md),
[`sbc()`](https://pedroliman.github.io/neuralsbi/reference/sbc.md),
[`tarp()`](https://pedroliman.github.io/neuralsbi/reference/tarp.md),
[`posterior_predictive()`](https://pedroliman.github.io/neuralsbi/reference/posterior_predictive.md)
– runs it through one execution path. That path is sequential by default
and parallel as soon as you declare a future plan:

    library(future)
    plan(multisession)          # or multicore, cluster, batchtools, ...
    fit <- npe(prior, simulator, n_simulations = 10000)

Nothing else changes: no argument to set, no variant of the function to
call. With no plan (or `plan(sequential)`) the simulations run in the
current process, and `neuralsbi` mentions the two lines above once per
session. Silence that hint with
`options(neuralsbi.parallel_hint = FALSE)`.

## Chunking

The parameter matrix is split into row chunks and one chunk at a time is
handed to the simulator, so a vectorized simulator keeps most of its
vectorization while progress can still be reported and work spread
across workers. The number of chunks depends only on the number of
simulations (`ceiling(n / chunk_size)`, with `chunk_size` chosen to give
about 64 chunks), *not* on the number of workers – which is what makes
results reproducible across plans. Set it yourself with the `chunk_size`
argument, or globally with `options(neuralsbi.chunks = )`. Larger chunks
mean less overhead per call and a coarser progress bar; smaller chunks
balance uneven simulator run times better. Under a multi-process plan
each chunk ships the simulator – and everything its environment captures
– to a worker, so a simulator closing over a large object is a reason to
raise `chunk_size`.

## Random numbers

Each chunk gets its own L'Ecuyer-CMRG random-number stream, derived from
the session's RNG state at the moment simulation starts. So
`set.seed(42)` before a fit (or `npe(..., seed = 42)`) gives the same
simulations whether you run sequentially or on 32 workers, and whatever
the worker count. What does change results is the chunking: fixing
`chunk_size` fixes the draws.

Because the simulations are drawn from separate streams, the simulator
no longer consumes the caller's RNG state directly; the state advances
by one draw per simulation phase.

## What is not parallelized

Training runs in the calling process. Torch tensors and modules are
external pointers that cannot be shipped to a worker, and libtorch
already uses multiple threads internally
([`torch::torch_set_num_threads()`](https://torch.mlverse.org/docs/reference/threads.html)).
The same goes for the posterior-sampling loops in
[`sbc()`](https://pedroliman.github.io/neuralsbi/reference/sbc.md) and
[`tarp()`](https://pedroliman.github.io/neuralsbi/reference/tarp.md) –
those call the trained network, not the simulator, and are reported with
their own progress bar.

## See also

[nsbi_progress](https://pedroliman.github.io/neuralsbi/reference/nsbi_progress.md)
for controlling progress bars.
