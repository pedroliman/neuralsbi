# Validate one distribution parameter of a named prior family

The constructors in
[prior_families](https://neuralsbi.pedrodelima.com/reference/prior_families.md)
take numeric vectors, one entry per parameter or one shared entry, and
every one of them ends up inside a `stats::d*()`/`q*()` pair. Those are
quiet about nonsense: `dgamma(x, shape = -1)` returns `NaN` with a
warning, `qgamma(u, shape = 0)` returns zeros, and either way the
mistake surfaces as a prior sample full of `NaN`s once the simulation
budget is already spent. The values are listed rather than described, as
[`check_counts()`](https://neuralsbi.pedrodelima.com/reference/check_counts.md)
does, since which entry is wrong is the thing worth reading.

## Usage

``` r
check_family_param(value, arg, positive = FALSE)
```

## Arguments

- value:

  The user's value.

- arg:

  Name of the argument.

- positive:

  Require every entry to be strictly positive, for a scale, a rate or a
  shape.

## Value

`value`, unchanged.

## Details

Names are left on the value:
[`family_param_names()`](https://neuralsbi.pedrodelima.com/reference/family_param_names.md)
reads the parameter names off whichever argument carries them, and that
has to happen after validation.
