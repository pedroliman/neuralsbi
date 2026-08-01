# Record one trace of the summed density at the shape of `tt`

Tracing refuses a captured tensor that carries a gradient, and the
fitted weights all do, so they are switched off around the recording and
switched back afterwards. Nothing else about the network changes.

## Usage

``` r
mdn_trace(de, xt, tt)
```
