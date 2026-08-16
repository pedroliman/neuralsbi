# `"resample"` starting points: pool prior draws, then SIR without replacement

A single pool of `max(n_pool, n_chains)` prior draws can land fewer than
`n_chains` of them in the posterior's support – the same shape of
problem
[`mcmc_init_proposal()`](https://neuralsbi.pedrodelima.com/reference/mcmc_init_proposal.md)
fixed for `"proposal"`. Padding the shortfall by recycling
([`rep_len()`](https://rdrr.io/r/base/rep.html) on too few indices) used
to start several chains from identical points, which defeats the
diagnostics downstream: split-Rhat and bulk ESS both assume the chains
they compare started from distinct locations, and a repeated start makes
chains look more agreed than they are.

## Usage

``` r
mcmc_init_resample(prior, log_prob_fn, n_chains, n_pool = 1000L)
```

## Details

Instead, more pools are drawn and their finite draws accumulated until
`n_chains` distinct ones are collected (up to 20 attempts, the same
budget
[`mcmc_init_proposal()`](https://neuralsbi.pedrodelima.com/reference/mcmc_init_proposal.md)
uses), and the final start is a weighted resample without replacement
from that accumulated, still-distinct set – not from a single pool. If
the budget runs out first, the error distinguishes a genuinely
degenerate surrogate (no draw at all landed in support) from a merely
low acceptance rate (some did, just not `n_chains` of them).
