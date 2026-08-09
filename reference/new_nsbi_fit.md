# Assemble the object an [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md), [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md) or [`nre()`](https://neuralsbi.pedrodelima.com/reference/nre.md) call returns

The three differ in the estimator they fit, the class they carry and a
field or two naming the architecture. Everything else – the fitted
estimator, the prior, the two standardizers, the dimensions, the names,
the simulation counts, the device – is the same list in the same order,
and every consumer downstream
([`save_npe()`](https://neuralsbi.pedrodelima.com/reference/save_npe.md),
[`check_fit_alive()`](https://neuralsbi.pedrodelima.com/reference/check_fit_alive.md),
[`cat_fit_common()`](https://neuralsbi.pedrodelima.com/reference/cat_fit_common.md),
[`summary()`](https://rdrr.io/r/base/summary.html)) reads it by name
from all three.

## Usage

``` r
new_nsbi_fit(de, prior, prep, class, extra = list())
```

## Arguments

- de:

  The fitted estimator. A ratio estimator travels in the same slot as a
  density estimator, since nothing downstream cares which it is until it
  comes time to score a `(theta, x)` pair.

- prior:

  The prior the simulations were drawn from.

- prep:

  The list from
  [`prepare_simulations()`](https://neuralsbi.pedrodelima.com/reference/prepare_simulations.md).

- class:

  The class to stamp on the result.

- extra:

  Named list of fields specific to this fit type, inserted before
  `device`.
