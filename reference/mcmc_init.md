# Starting points for the chains

`"resample"` (the default, and `sbi`'s) draws a large pool from the
prior, weights it by the posterior density and resamples without
replacement: a sampling-importance-resampling start that puts the chains
where the mass is, which matters because slice sampling has no
adaptation phase to rescue a bad start. `"proposal"` just takes prior
draws.

## Usage

``` r
mcmc_init(
  prior,
  log_prob_fn,
  n_chains,
  strategy = c("resample", "proposal"),
  n_pool = 1000L
)
```

## Details

Two departures from `sbi`, both about the case this exists for. `sbi`
draws its pool once per chain, so 20 chains cost 20 pools; one pool
shared across chains costs a twentieth of that and, drawn without
replacement, gives distinct starting points rather than independent ones
that may coincide. And the weighting stays in log space: `sbi`
normalizes the log weights and exponentiates, which is fine until a few
thousand independent observations spread the log-likelihood over prior
draws by thousands of nats and every weight but a handful underflows to
zero.
