# Validate a vector of probabilities

[`check_prob()`](https://neuralsbi.pedrodelima.com/reference/check_prob.md)
for an argument that is a vector rather than one number, such as the
credibility `levels` of
[`expected_coverage()`](https://neuralsbi.pedrodelima.com/reference/expected_coverage.md).
A level outside `(0, 1)` produces one meaningless number in a table of
plausible ones: `levels = c(-1, 2)` scores the empty interval and the
whole line, and comes back as coverage 0 and 1 with nothing said. The
message lists the values, as
[`check_counts()`](https://neuralsbi.pedrodelima.com/reference/check_counts.md)
does, since which entry is wrong is the thing worth reading.

## Usage

``` r
check_probs(p, arg, open = TRUE)
```

## Arguments

- p:

  The user's value.

- arg:

  Name of the argument.

- open:

  Require `0 < p < 1` (the default). `FALSE` allows the endpoints.

## Value

`p` as a double vector.
