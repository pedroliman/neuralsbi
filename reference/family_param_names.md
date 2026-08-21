# Parameter names taken from whichever argument carries them

[`prior_uniform()`](https://neuralsbi.pedrodelima.com/reference/prior_uniform.md)
reads them off `low` or `high`; a family with three arguments has three
places to look, and the first complete set wins.

## Usage

``` r
family_param_names(args, d)
```

## Arguments

- args:

  Named list of the raw (un-recycled) arguments.

- d:

  Number of parameters.
