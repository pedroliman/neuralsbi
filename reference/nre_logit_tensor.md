# Per-row logit of the classifier, as a torch tensor

This is the log ratio itself: the module's single output for each
`(theta, x)` row.

## Usage

``` r
nre_logit_tensor(net, theta, x)
```
