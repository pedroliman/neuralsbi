# Describe one parameter set in an error message

The values that produced a failure, as `mu = 1.83, sigma = 0.42`, or
bare numbers when the prior does not name its parameters. Rounded,
because four significant digits are enough to recognise a draw and a
full-precision double is not. Only the first `max_show` are printed so a
40-parameter model does not fill the console.

## Usage

``` r
describe_params(theta_i, max_show = 6L)
```

## Arguments

- theta_i:

  One parameter set, named or not.

- max_show:

  How many values to print before truncating.
