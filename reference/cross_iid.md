# Walk the theta blocks of the cross product, calling `collect(idx, lp)`

`lp` is the flat vector `score()` returns for the block, running
theta-major: the whole observation set for the first parameter, then for
the second, and so on. `score` is
[de_log_prob()](https://neuralsbi.pedrodelima.com/reference/density_estimator.md)
for a density estimator and
[`nre_score()`](https://neuralsbi.pedrodelima.com/reference/nre_score.md)
for a ratio estimator; the blocking is the same either way.

## Usage

``` r
cross_iid(de, x, theta, max_batch, collect, score = de_log_prob)
```

## Details

`theta` is blocked first, exactly as before. What changed (#248) is that
a theta block used to hand `score()` its whole observation set in one
call – fine when `n_obs <= max_batch`, but when `n_obs` alone exceeds
`max_batch` (NLE/NRE's headline case of conditioning on thousands of
trials) that one call saw `n_obs` pairs regardless of `max_batch`, the
observation side never being chunked at all. Observations are now
chunked within each theta block too, so no call to `score()` sees more
than `max_batch` pairs.
[`collect()`](https://dplyr.tidyverse.org/reference/compute.html) still
only sees one call per theta block, with the same complete
`length(idx) * n_obs` flat vector it always did – the inner chunking is
accumulated into that vector before
[`collect()`](https://dplyr.tidyverse.org/reference/compute.html) runs,
so callers are untouched.
[`mdn_iid_blocks()`](https://neuralsbi.pedrodelima.com/reference/mdn_iid_blocks.md)
already chunked both sides for the MDN path; this brings the flow/NRE
path in line with it.
