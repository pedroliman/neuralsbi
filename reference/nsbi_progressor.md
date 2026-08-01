# Create a progress reporter for one phase

Returns `function(amount = 1, total = NULL, done = FALSE)`. `total`
revises the denominator mid-flight, which is what training needs: the
epoch a run will stop at is only known once it stops. Under progressr –
whose step count is fixed at creation – the revised fraction is mapped
back onto the original `steps`, so both reporters show the same picture.

## Usage

``` r
nsbi_progressor(steps, label = NULL)
```

## Arguments

- steps:

  Expected number of units of work.

- label:

  Short phase name shown next to the bar.
