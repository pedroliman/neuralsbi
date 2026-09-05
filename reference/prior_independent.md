# Combine independent priors into one joint prior

Stan writes a joint prior as one sampling statement per parameter and
lets the product take care of itself. This is that, as an object: give
it the per-parameter priors in the order the simulator expects them and
it returns the product prior, with `dim` the total, `lower`/`upper`
stacked, and a `log_prob` that sums the parts. It is the thing most
[`prior_custom()`](https://neuralsbi.pedrodelima.com/reference/prior_custom.md)
calls were written to do by hand.

## Usage

``` r
prior_independent(...)
```

## Arguments

- ...:

  Priors, one per block of parameters. Name them
  (`prior_independent(beta = ..., gamma = ...)`) to name the parameters:
  a name given here beats the component's own `param_names` for a
  one-parameter component, and a multi-parameter component keeps its own
  names, or takes `name1`, `name2`, ... when it has none.

## Value

An `nsbi_prior` object.

## Details

Components may themselves cover several parameters, so
`prior_independent(prior_normal(mean = c(0, 0)), prior_beta(2, 15))` is
a three-parameter prior. Any `nsbi_prior` works as a component,
including a
[`prior_custom()`](https://neuralsbi.pedrodelima.com/reference/prior_custom.md);
the result can only be written out by
[`stan_code()`](https://neuralsbi.pedrodelima.com/reference/stan_export.md)
when every component comes from a named family (see
[prior_families](https://neuralsbi.pedrodelima.com/reference/prior_families.md)).

## See also

[prior_families](https://neuralsbi.pedrodelima.com/reference/prior_families.md)
for the components,
[`prior_truncated()`](https://neuralsbi.pedrodelima.com/reference/prior_truncated.md)
to bound one,
[priors](https://neuralsbi.pedrodelima.com/reference/priors.md) for the
rest.

## Examples

``` r
prior <- prior_independent(
  p_S1S2 = prior_beta(2, 15),
  hr_S1  = prior_lognormal(log(3), 0.3),
  hr_S2  = prior_lognormal(log(10), 0.25)
)
prior
#> <nsbi_prior> type=independent, dim=3
#>   parameters: p_S1S2, hr_S1, hr_S2 
#>   lower: 0, 0, 0 
#>   upper: 1, Inf, Inf 
#>   p_S1S2 ~ beta(2, 15)
#>   hr_S1 ~ lognormal(1.099, 0.3)
#>   hr_S2 ~ lognormal(2.303, 0.25)
sample_prior(prior, 3)
#>         p_S1S2    hr_S1    hr_S2
#> [1,] 0.1385162 2.483443 12.77441
#> [2,] 0.1078683 1.788164  6.66167
#> [3,] 0.1462143 2.657680 10.04431
```
