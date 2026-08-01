# Is this torch module still backed by a live pointer?

[`readRDS()`](https://rdrr.io/r/base/readRDS.html) on a torch-backed fit
returns a module whose external pointer is nil. Nothing about the R
object says so, so the only way to find out is to touch a tensor.

## Usage

``` r
torch_net_alive(net)
```
