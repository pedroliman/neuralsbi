# Replay the MDN's i.i.d. density as TorchScript instead of driving it from R

Every operation in
[`mdn_iid_blocks()`](https://neuralsbi.pedrodelima.com/reference/mdn_iid_blocks.md)
crosses from R into libtorch, and at MCMC batch sizes that crossing
costs far more than the arithmetic behind it: about 0.2 ms each, thirty
of them per evaluation, against a few hundred microseconds of actual
work.
[`torch::jit_trace()`](https://torch.mlverse.org/docs/reference/jit_trace.html)
records the same computation once and replays it inside libtorch, so an
evaluation costs one crossing rather than thirty.

## Usage

``` r
mdn_trace_cache(de, xt, max_batch, eager, warmup = 4L)
```

## Arguments

- de, xt, max_batch:

  As in
  [`mdn_iid_blocks()`](https://neuralsbi.pedrodelima.com/reference/mdn_iid_blocks.md),
  with the observations already a tensor.

- eager:

  The evaluator to check each trace against, and to fall back to.

- warmup:

  Calls to serve eagerly before recording anything.

## Value

`function(theta)` returning a traced function for that many rows, or
`NULL` when the eager path should be used. `NULL` if tracing is switched
off.

## Details

It is the same code either way. Tracing runs the eager path and records
what it did, so there is no separate implementation to keep in step with
the eager one – that claim holds only for the traced path relative to
[`mdn_log_prob_tensor()`](https://neuralsbi.pedrodelima.com/reference/mdn_log_prob_tensor.md).
The MDN density still has three implementations in total, one per
runtime:
[`mdn_log_prob_tensor()`](https://neuralsbi.pedrodelima.com/reference/mdn_log_prob_tensor.md)
(the eager training path),
[`mdn_mixture()`](https://neuralsbi.pedrodelima.com/reference/mdn_mixture.md)/[`mdn_chunk_lp()`](https://neuralsbi.pedrodelima.com/reference/mdn_chunk_lp.md)
(the i.i.d. fast path used here), and `stan_fn_mdn()` (the generated
Stan code, `R/stan.R`). That is by design, not drift: three runtimes
need three implementations, and the tests pin them to each other
numerically.

Three things make this a shortcut rather than the path. Recording a
trace costs several evaluations' worth of time, so nothing is recorded
until the evaluator has been called `warmup` times and it is clear this
is a loop rather than a one-off; a single
[`log_lik()`](https://neuralsbi.pedrodelima.com/reference/log_lik.md)
call should not pay for a compiler. A trace fixes the shapes it was
recorded at, so there is one per number of parameter rows, checked
against the eager result before anything uses it. And it is only worth
recording when the whole observation set fits in one chunk, since
otherwise the graph unrolls the chunk loop. Failing any of these is not
an error: the caller falls back to the eager path, which is why `NULL`
is a perfectly good answer here.

Set `options(neuralsbi.jit = FALSE)` to skip tracing entirely.
