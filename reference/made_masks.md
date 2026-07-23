# MADE masks for one autoregressive transform.

Degrees: theta inputs get 1..p, conditioning inputs get 0 (visible to
all), hidden units cycle through 1..(p-1) (or 0 when p = 1). A
connection into a hidden unit requires hidden_degree \>= input_degree; a
connection into output dimension d requires d \> hidden_degree. This
makes output d a function of \\\theta\_{\<d}\\ and x only.

## Usage

``` r
made_masks(dim_theta, dim_x, hidden)
```
