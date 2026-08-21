# Priors for neural simulation-based inference

A prior in `neuralsbi` is a lightweight object (class `nsbi_prior`) that
knows how to (a) draw samples and (b) evaluate its log-density. Bounded
priors also carry `lower`/`upper` support limits, which are used to
reject out-of-support posterior samples ("leakage" correction).

## Details

There are four ways to build one.
[`prior_uniform()`](https://neuralsbi.pedrodelima.com/reference/prior_uniform.md)
is the box prior most benchmark tasks use.
[`prior_normal()`](https://neuralsbi.pedrodelima.com/reference/prior_normal.md)
and the named families in
[prior_families](https://neuralsbi.pedrodelima.com/reference/prior_families.md)
(log-normal, exponential, gamma, beta, Student-t, Cauchy and the half
versions of the last two) give one marginal per parameter, under Stan's
argument names.
[`prior_independent()`](https://neuralsbi.pedrodelima.com/reference/prior_independent.md)
multiplies those together into a joint prior, and
[`prior_truncated()`](https://neuralsbi.pedrodelima.com/reference/prior_truncated.md)
bounds one, renormalizing the density by the mass it keeps.
[`prior_custom()`](https://neuralsbi.pedrodelima.com/reference/prior_custom.md)
takes a sampler and a density you write yourself, for anything the
families do not cover.

Prefer a named family over
[`prior_custom()`](https://neuralsbi.pedrodelima.com/reference/prior_custom.md)
where one fits. A family carries its support bounds into the leakage
correction, survives
[`prior_truncated()`](https://neuralsbi.pedrodelima.com/reference/prior_truncated.md)
with an exact normalizing constant, and is what
[`stan_code()`](https://neuralsbi.pedrodelima.com/reference/stan_export.md)
writes out as a sampling statement; a custom prior does none of those.

## See also

[prior_families](https://neuralsbi.pedrodelima.com/reference/prior_families.md),
[`prior_independent()`](https://neuralsbi.pedrodelima.com/reference/prior_independent.md),
[`prior_truncated()`](https://neuralsbi.pedrodelima.com/reference/prior_truncated.md),
[`sample_prior()`](https://neuralsbi.pedrodelima.com/reference/sample_prior.md),
[`within_support()`](https://neuralsbi.pedrodelima.com/reference/within_support.md).
