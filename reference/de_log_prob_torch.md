# Shared tensor plumbing behind every neural `de_log_prob.*` method

Coerces `theta` and `x` to matrices, broadcasts a single-row `x` up to
`theta`'s row count (the same broadcast
[`lingauss_mean()`](https://neuralsbi.pedrodelima.com/reference/lingauss_mean.md)'s
caller does on `mu`, just on the other operand), moves both to torch,
and evaluates `log_prob_fn` under `with_no_grad()`. `log_prob_fn` is the
per-estimator tensor function –
[`mdn_log_prob_tensor()`](https://neuralsbi.pedrodelima.com/reference/mdn_log_prob_tensor.md),
[`maf_log_prob_tensor()`](https://neuralsbi.pedrodelima.com/reference/maf_log_prob_tensor.md)
or
[`nsf_log_prob_tensor()`](https://neuralsbi.pedrodelima.com/reference/nsf_log_prob_tensor.md).

## Usage

``` r
de_log_prob_torch(de, theta, x, log_prob_fn)
```
