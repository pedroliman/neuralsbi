# How this simulator receives its parameters

Decided once per run from
[`formals()`](https://rdrr.io/r/base/formals.html), never by probing the
simulator: a wrong guess would produce a wrong posterior with no error.

## Usage

``` r
sim_dispatch(simulator, param_names)
```

## Value

`"named"` (one scalar per formal) or `"vector"` (the named parameter
vector as the first argument).

## Details

A partial match warns. When some parameter names appear among the
formals and others do not, both signatures are plausible and the vector
form wins, which sends the whole parameter vector to the first formal
and lets the rest take their defaults. That trains on nonsense without
erroring, and one typo in a prior name is enough to cause it. A warning
rather than an error, because a vector-signature simulator whose first
argument happens to carry a parameter's name is legal.
