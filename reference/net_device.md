# The torch device a fitted net's parameters currently live on

`de_log_prob()`/`de_sample()` need to know where to put their input
tensors, and the one place that is always true is the net itself –
`de$device` (set at fit time) is a record of that, not a guarantee,
since nothing stops a user from moving the net with `net$to()`
afterwards. A net with no parameters (there is no such estimator today,
but nothing rules one out) defaults to CPU.

## Usage

``` r
net_device(net)
```
