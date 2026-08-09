# Per-column scale a relative ridge is measured against

One scale per column, so a matrix whose columns span several orders of
magnitude is regularized by the same *relative* amount everywhere. A
single shared scale – the mean or the maximum of the diagonal – is no
better than an absolute ridge here: it is set by the largest column and
swamps the smallest.

## Usage

``` r
ridge_scale(scale)
```

## Arguments

- scale:

  Per-column scales, in the units of the matrix diagonal.

## Details

A column with no scale of its own borrows the largest one that has any,
and if nothing does the scale is 1 and the ridge acts absolutely. That
is the degenerate case the ridge exists for, and it is the only one
where the answer is arbitrary.
