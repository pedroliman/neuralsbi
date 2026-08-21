# Assemble marginals into an independent joint prior

Assemble marginals into an independent joint prior

## Usage

``` r
marginal_prior(marginals, param_names = NULL, type = "independent")
```

## Arguments

- marginals:

  List of marginals, one per parameter.

- param_names:

  Optional parameter names.

- type:

  The `type` recorded on the prior.

## Value

An `nsbi_prior` object.
