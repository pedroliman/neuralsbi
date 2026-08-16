# `"proposal"` starting points: pool prior draws until enough land in support

A single batch of `n_chains` prior draws is finite only with probability
`acceptance^n_chains`, so a posterior that is entirely ordinary but
excludes a third of the prior's support (acceptance 0.65) already fails
almost every batch: `0.65^20` is about 3e-4. Requiring one batch to
succeed whole conflates that low-but-fine acceptance rate with true
degeneracy. Drawing a pool per attempt and keeping every finite draw
across attempts, the way `"resample"` already pools, fixes this: the
chance of collecting `n_chains` finite draws now grows with the *total*
number sampled, not with the size of one lucky batch.

## Usage

``` r
mcmc_init_proposal(prior, log_prob_fn, n_chains, n_pool = 1000L)
```

## Details

Up to 20 pools of `max(n_pool, n_chains)` prior draws are drawn
(`n_pool` defaults to 1000, so 20000 draws by default), stopping as soon
as enough finite ones have accumulated. If the budget runs out first,
the error distinguishes two cases: no draw at all landed in support,
which is consistent with a degenerate surrogate, versus some draws did
but not enough, which is just a low acceptance rate and should not be
reported as degeneracy.
