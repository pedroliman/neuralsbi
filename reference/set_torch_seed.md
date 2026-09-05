# Seed torch's global RNG, returning its previous state to restore later

[`torch::torch_manual_seed()`](https://torch.mlverse.org/docs/reference/torch_manual_seed.html)
reseeds torch's one global generator, and there is no way to ask it for
an independent stream instead. A caller that reseeds it and never
restores the previous state leaves every later *unseeded* torch call in
the same session – another fit, another `de_sample()` call, a
[`c2st()`](https://neuralsbi.pedrodelima.com/reference/c2st.md) call
with no `seed` of its own – drawing from wherever the seeded call left
the generator, rather than from a fresh one (GitHub \#275, the torch
analogue of \#272).

## Usage

``` r
set_torch_seed(seed)
```

## Arguments

- seed:

  Integer seed to set torch's global RNG to.

## Value

The RNG state captured before reseeding, suitable for
[`torch::torch_set_rng_state()`](https://torch.mlverse.org/docs/reference/torch_manual_seed.html).

## Details

Pair this with `on.exit(torch::torch_set_rng_state(old), add = TRUE)`
right after calling it. That is the same save-then-restore shape
[`with_fixed_seed()`](https://neuralsbi.pedrodelima.com/reference/with_fixed_seed.md)
uses for R's own RNG, but split into a value-returning helper plus the
caller's own [`on.exit()`](https://rdrr.io/r/base/on.exit.html) rather
than an expression-wrapping one: both call sites need the seed live for
the rest of a long-running function – a network's random initialization
and the epochs of shuffled batches that follow it, not just one
expression's value – so the restore has to fire when that function
returns, not when this helper does.

Requires `torch` (call
[`require_torch()`](https://neuralsbi.pedrodelima.com/reference/require_torch.md)
first; both call sites are already past that check when they reach
this).
