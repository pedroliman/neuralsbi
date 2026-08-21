# Stan's `T[,]` clause for a marginal cut down from its natural support

Written out even though fixed bounds and fixed parameters make it a
constant that Stan would drop. The generated program is meant to be read
and edited, and a reader who lifts the prior parameters into `data`
would otherwise inherit a silently wrong normalization.

## Usage

``` r
stan_truncation(m)
```

## Arguments

- m:

  One marginal.
