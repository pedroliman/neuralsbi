# Seed R's RNG for the rest of the calling function, then restore it

Entry points with a `seed` argument call this once near the top. It does
`set.seed(seed)` and registers an
[`on.exit()`](https://rdrr.io/r/base/on.exit.html) in the *caller's*
frame that puts the caller's `.Random.seed` back (or removes it if there
was none), so a seeded call is reproducible without leaving the user's
own random stream in a state fixed by `seed` (GitHub \#381, the
entry-point analogue of \#272 and \#319). With `seed = NULL` it does
nothing. Use
[`with_fixed_seed()`](https://neuralsbi.pedrodelima.com/reference/with_fixed_seed.md)
instead when only one expression needs the seed.

## Usage

``` r
local_seed(seed, envir = parent.frame())
```

## Arguments

- seed:

  Integer seed, or `NULL` to leave the RNG alone.

- envir:

  Frame whose exit triggers the restore.

## Value

`NULL`, invisibly.
