# Draw posterior samples for one diagnostic trial, insisting on the full count

[`sample.nsbi_posterior()`](https://neuralsbi.pedrodelima.com/reference/sample.nsbi_posterior.md)
returns fewer rows than asked for when a bounded prior and a leaky
estimator leave rejection sampling short, and only warns. The
diagnostics cannot absorb that quietly.
[`sbc()`](https://neuralsbi.pedrodelima.com/reference/sbc.md) bins its
ranks against `n_posterior_samples`, so a trial that came back short is
scored on a scale it was never drawn on, and the ranks compress toward
zero; rescaling that one trial on its own would instead make it
incomparable to the others. Either way the run reports a miscalibrated
posterior when the real cause is lost draws, so stop and say so.

## Usage

``` r
diagnostic_draws(post, n, trial)
```
