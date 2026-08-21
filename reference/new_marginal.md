# One parameter's marginal: a family, its parameters, and its bounds

The canonical form every named-family prior is stored in, under
`prior$params$marginals`. `lower`/`upper` are the effective bounds, that
is the family's own support intersected with whatever truncation was
asked for, and `log_norm` is the log of the probability mass left inside
them. Holding the truncated mass here rather than recomputing it per
call is what keeps the density proper without paying for a CDF
evaluation on every row.

## Usage

``` r
new_marginal(family, args, lower = -Inf, upper = Inf, label = "the parameter")
```

## Arguments

- family:

  Family name, a key of
  [`prior_family()`](https://neuralsbi.pedrodelima.com/reference/prior_family.md).

- args:

  Named list of scalar distribution parameters.

- lower, upper:

  Truncation bounds, before intersecting with the family's own support.

- label:

  How to refer to this parameter in an error message.

## Value

A marginal specification.
