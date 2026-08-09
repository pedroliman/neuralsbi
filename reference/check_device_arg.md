# Validate the `device` argument shared by [`npe()`](https://neuralsbi.pedrodelima.com/reference/npe.md) and [`nle()`](https://neuralsbi.pedrodelima.com/reference/nle.md)

Only checks that `device` is one of the recognized keywords – `"cpu"`,
`"cuda"`, `"mps"`, or `"gpu"`/`"auto"` – so a typo is caught before the
simulator runs, the same reason
[`check_train_controls()`](https://neuralsbi.pedrodelima.com/reference/check_train_controls.md)
runs early. This is the whole check for `"cpu"`. Turning
`"cuda"`/`"mps"`/`"gpu"`/`"auto"` into a concrete, available device
needs `torch` loaded, which `density_estimator = "linear_gaussian"`
never requires; that step is
[`resolve_device()`](https://neuralsbi.pedrodelima.com/reference/resolve_device.md)'s
job, called only once a `torch`-backed estimator is about to train.

## Usage

``` r
check_device_arg(device)
```

## Arguments

- device:

  The user's value.

## Value

`device`, unchanged.
