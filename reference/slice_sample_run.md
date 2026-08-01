# The loop behind [`slice_sample()`](https://neuralsbi.pedrodelima.com/reference/slice_sample.md), split out so the progress context wraps it

The loop behind
[`slice_sample()`](https://neuralsbi.pedrodelima.com/reference/slice_sample.md),
split out so the progress context wraps it

## Usage

``` r
slice_sample_run(
  log_prob_fn,
  init,
  n_draws,
  warmup,
  thin,
  width,
  max_steps,
  verbose
)
```
