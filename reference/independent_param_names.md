# Parameter names for a product prior

An argument name wins for a one-parameter component, since that is the
obvious way to write it; a wider component keeps the names it already
has, and otherwise the argument name is numbered. Anything short of a
complete set of names is dropped, because a partly named prior would
label some columns and not others.

## Usage

``` r
independent_param_names(parts, dims, nms)
```

## Arguments

- parts:

  The component priors.

- dims:

  Their dimensions.

- nms:

  The names of the `...` arguments, or `NULL`.
