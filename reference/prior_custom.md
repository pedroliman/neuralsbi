# Build a prior from arbitrary sampling / density functions

This is the one prior a user writes by hand, so it is the one where a
mistake is most likely, and the mistakes are quiet. A `log_prob_fn` that
returns a single number instead of one per row is legal R that surfaces
much later as an MCMC initialization failure; a `lower` of the wrong
length is recycled by [`sweep()`](https://rdrr.io/r/base/sweep.html)
into a support test that rejects the wrong draws. Everything is
therefore checked at construction, including one probe call of
`sample_fn(2)` and, when it is given, `log_prob_fn()` on those two rows.

## Usage

``` r
prior_custom(
  sample_fn,
  log_prob_fn = NULL,
  dim,
  lower = NULL,
  upper = NULL,
  param_names = NULL
)
```

## Arguments

- sample_fn:

  Function `function(n)` returning an `n x dim` matrix, one row per
  draw.

- log_prob_fn:

  Function `function(theta)` returning a length-`n` vector of log
  densities, one per row of `theta`. Optional; required only for
  methods/diagnostics that need it.

- dim:

  Number of parameters.

- lower, upper:

  Optional support bounds enabling out-of-support rejection. Numeric of
  length `dim`, or length 1 to apply the same bound to every parameter.

- param_names:

  Optional character vector of parameter names, one per parameter. These
  name the columns of every downstream parameter matrix, posterior
  sample and diagnostic plot, and they decide how the simulator is
  called: a simulator whose formals are exactly these names receives one
  scalar per formal, and any other simulator receives the parameter
  vector as its first argument. See
  [nsbi_simulator](https://neuralsbi.pedrodelima.com/reference/nsbi_simulator.md).

## Value

An `nsbi_prior` object.

## Stan export

A custom prior is arbitrary R code rather than a named distribution with
parameters, so
[`stan_code()`](https://neuralsbi.pedrodelima.com/reference/stan_export.md)
cannot restate it as a Stan sampling statement. Take
`stan_code(fit, model = FALSE)` and write the model block yourself, or
build the prior out of
[prior_families](https://neuralsbi.pedrodelima.com/reference/prior_families.md)
and
[`prior_independent()`](https://neuralsbi.pedrodelima.com/reference/prior_independent.md)
instead, which
[`stan_code()`](https://neuralsbi.pedrodelima.com/reference/stan_export.md)
does write out.

## Examples

``` r
prior <- prior_custom(
  sample_fn = function(n) cbind(rexp(n, 1), rexp(n, 2)),
  log_prob_fn = function(theta) {
    dexp(theta[, 1], 1, log = TRUE) + dexp(theta[, 2], 2, log = TRUE)
  },
  dim = 2, lower = 0, param_names = c("beta", "gamma")
)
```
