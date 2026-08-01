# Pick a progressr handler that can show an ETA

`handler_txtprogressbar` is the progressr default and reports no ETA, so
prefer the cli handler when cli is available. A user who has registered
handlers of their own never gets here.

## Usage

``` r
default_progressr_handler()
```
