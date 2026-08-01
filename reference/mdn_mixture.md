# The mixture an MDN puts over its target for each row of `tt`

Split out from
[`mdn_iid_blocks()`](https://neuralsbi.pedrodelima.com/reference/mdn_iid_blocks.md)
so the eager driver and the traced one share a single definition of the
density. Everything here depends on the conditioning variable alone, so
it is computed once per call and reused across every observation chunk.

## Usage

``` r
mdn_mixture(de, tt)
```
