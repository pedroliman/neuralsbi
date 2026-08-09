# Resolve a validated `device` keyword to a concrete, available torch device

[`check_device_arg()`](https://neuralsbi.pedrodelima.com/reference/check_device_arg.md)
has already rejected anything but `"cpu"`, `"cuda"`, `"mps"`, `"gpu"` or
`"auto"`; this turns that keyword into `"cpu"`, `"cuda"` or `"mps"`,
matching what Python `sbi`'s `process_device()` does. Requires `torch`
(call
[`require_torch()`](https://neuralsbi.pedrodelima.com/reference/require_torch.md)
first; every caller here is about to build a network, so it already
has).

## Usage

``` r
resolve_device(device)
```

## Arguments

- device:

  One of `"cpu"`, `"cuda"`, `"mps"`, `"gpu"`, `"auto"`.

## Value

`"cpu"`, `"cuda"` or `"mps"`.

## Details

`"cuda"`/`"mps"` name a specific device, so asking for one that is not
there errors rather than downgrading silently – a silent fallback would
hide a real problem (a missing CUDA build, a non-Apple-silicon Mac).
`"gpu"`/`"auto"` never named a specific device, so it resolves CUDA -\>
MPS -\> CPU and *can* fall back silently, mirroring `sbi`'s `"gpu"`.
