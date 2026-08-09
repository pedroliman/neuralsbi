# Shared tensor plumbing behind `de_sample.nsbi_de_maf` and `de_sample.nsbi_de_nsf`

Takes the single conditioning row, replicates it to `n` rows, draws a
standard-normal base sample, and inverts the flow – all on the net's
device (see
[`net_device()`](https://neuralsbi.pedrodelima.com/reference/net_device.md)),
then brings the draws back to CPU as a plain R matrix. `inverse_fn` is
the per-flow inverse –
[`maf_inverse()`](https://neuralsbi.pedrodelima.com/reference/maf_inverse.md)
or
[`nsf_inverse()`](https://neuralsbi.pedrodelima.com/reference/nsf_inverse.md).
The MDN has no inverse to share here; it samples its mixture directly.

## Usage

``` r
de_sample_flow(de, x, n, inverse_fn)
```
