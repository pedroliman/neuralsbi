# MCMC over a learned likelihood

An [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) fit
gives an unnormalized posterior, \\q\_\phi(x \mid \theta)\\p(\theta)\\,
but no way to draw from it directly. `neuralsbi` samples it with a
univariate slice sampler (Neal, 2003) run over many chains at once.

## Details

Slice sampling is the default in Python `sbi` for the same reason it is
the default here: it has no step size to tune, adapts its scale to the
target as it goes, and handles bounded supports without any special
casing, because a prior that returns `-Inf` outside its support simply
shrinks the slice interval.

The vectorization runs across chains rather than across dimensions.
Every chain proposes its next value for the same coordinate at the same
time, so one step costs one batched call to the density estimator
instead of `n_chains` separate forward passes. With a neural likelihood
that difference is the whole running time.

## References

Neal, R. M. (2003). Slice sampling. *The Annals of Statistics* 31(3),
705-767.

## See also

[`posterior.nsbi_nle()`](https://neuralsbi.pedrodelima.com/reference/posterior.nsbi_nle.md)
for the arguments that control it, and
[`stan_code()`](https://neuralsbi.pedrodelima.com/reference/stan_export.md)
for handing the same likelihood to NUTS instead.
