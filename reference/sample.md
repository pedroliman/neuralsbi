# Draw samples (S3 generic)

`neuralsbi` turns [`base::sample()`](https://rdrr.io/r/base/sample.html)
into an S3 generic so that `sample(posterior, n)` reads the way
statisticians expect. For any object without a dedicated method
(vectors, etc.) this falls back to
[`base::sample()`](https://rdrr.io/r/base/sample.html) unchanged.

## Usage

``` r
sample(x, ...)

# Default S3 method
sample(x, ...)
```

## Arguments

- x:

  Object to sample from.

- ...:

  Passed on to methods /
  [`base::sample()`](https://rdrr.io/r/base/sample.html).
