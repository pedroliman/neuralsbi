# Apply `fun` to each batch, sequentially or across future workers

The parallel branch keeps at most one future per worker in flight and
polls for completions, so progress is reported as batches finish rather
than in one jump at the end. `fun` is called as `fun(batch, tick)`:
running sequentially it ticks the bar itself, once per unit of work;
running on a worker it cannot report back mid-batch, so `tick` does
nothing and the batch is credited on completion.

## Usage

``` r
nsbi_batch_apply(batches, fun, p = NULL, weights = NULL, seeds = NULL)
```

## Arguments

- batches:

  List of arguments to `fun`.

- fun:

  Function of a batch and a `tick()` callback.

- p:

  Progress reporter from
  [`nsbi_progressor()`](https://neuralsbi.pedrodelima.com/reference/nsbi_progressor.md),
  or `NULL`.

- weights:

  Progress units to credit per batch (defaults to 1 each).

- seeds:

  L'Ecuyer-CMRG seeds, one per batch, handed to the future backend.
  `fun` is free to set its own streams inside the batch.
