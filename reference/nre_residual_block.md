# One pre-activation residual block, as in `nflows`' `ResidualNet`

`relu -> linear -> relu -> linear`, added back to the input. The second
linear layer starts near zero so the block starts as the identity, the
same trick
[`made_module()`](https://neuralsbi.pedrodelima.com/reference/made_module.md)
uses on its output heads and the reason a deep stack trains at all from
a cold start.

## Usage

``` r
nre_residual_block(features)
```
