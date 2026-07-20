# Build a prior from arbitrary sampling / density functions

Build a prior from arbitrary sampling / density functions

## Usage

``` r
prior_custom(sample_fn, log_prob_fn = NULL, dim, lower = NULL, upper = NULL)
```

## Arguments

- sample_fn:

  Function `function(n)` returning an `n x dim` matrix.

- log_prob_fn:

  Function `function(theta)` returning a length-`n` vector of log
  densities. Optional; required only for methods/diagnostics that need
  it.

- dim:

  Number of parameters.

- lower, upper:

  Optional support bounds (numeric vectors) enabling out-of-support
  rejection.

## Value

An `nsbi_prior` object.
