# Should *we* draw a bar?

Distinct from
[`progress_backend()`](https://neuralsbi.pedrodelima.com/reference/progress_backend.md):
progressr users who registered their own handlers always get updates,
whatever this returns. This governs only whether `neuralsbi` installs a
handler itself (or draws its own bar).

## Usage

``` r
progress_render()
```
