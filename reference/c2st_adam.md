# Adam for [`c2st_mlp_prob()`](https://neuralsbi.pedrodelima.com/reference/c2st_mlp_prob.md), preferring torch's C++ implementation

`optim_ignite_adam()` takes the whole step in C++ where `optim_adam()`
loops over the parameters in R. The two produce identical trajectories,
and this loop is thousands of steps over six small tensors, so the
R-side overhead is most of what it costs: swapping one for the other
halves the running time. `optim_ignite_adam()` arrived in torch 0.13.0
and DESCRIPTION allows 0.11.0, which is what the fallback is for.

## Usage

``` r
c2st_adam(groups, lr)
```

## Arguments

- groups:

  Parameter groups, as `optim_adam()` takes them.

- lr:

  Adam step size.
