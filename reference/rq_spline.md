# Monotonic rational-quadratic spline, batched.

Monotonic rational-quadratic spline, batched.

## Usage

``` r
rq_spline(
  inputs,
  w_un,
  h_un,
  d_un,
  inverse = FALSE,
  tail_bound = 3,
  min_bin = 0.001,
  min_deriv = 0.001
)
```

## Arguments

- inputs:

  `(N,)` tensor of values to transform.

- w_un, h_un, d_un:

  Unnormalized widths `(N, K)`, heights `(N, K)`, and interior
  derivatives `(N, K - 1)`.

- inverse:

  Apply the inverse transform.

- tail_bound:

  Spline acts on `[-tail_bound, tail_bound]`; identity outside.

## Value

`list(outputs, logdet)`, both `(N,)` tensors.
