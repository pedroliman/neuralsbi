# Progress reporting

Long-running steps – running the simulator, training a neural density
estimator – report progress through one mechanism, so a fit that
simulates and then trains shows the same style of bar for both phases.
Each bar reports the elapsed time and an ETA extrapolated from the work
done so far.

## Details

If the progressr package is installed, `neuralsbi` emits standard
progressr updates, so
[`progressr::handlers()`](https://progressr.futureverse.org/reference/handlers.html)
and
[`progressr::with_progress()`](https://progressr.futureverse.org/reference/with_progress.html)
control reporting exactly as they do for any other progressr-aware
package:

    library(progressr)
    handlers(global = TRUE)
    handlers("cli")            # or "progress", "txtprogressbar", "rstudio"
    fit <- npe(prior, simulator, n_simulations = 10000)

Without progressr installed (or when it is installed but no handler is
registered) `neuralsbi` draws its own bar, which needs no extra
packages.

Reporting is controlled by `options(neuralsbi.progress = )`:

- `"auto"` (default) – report in interactive sessions, stay quiet in
  scripts, `R CMD check`, and knitr.

- `TRUE` – always report, interactive or not.

- `FALSE` – never report.

- `"builtin"` – always report with the built-in bar, ignoring progressr.

## Training progress

While a neural estimator trains, one progress step is one epoch. The
total is not known in advance because training stops early, so the bar
targets the epoch at which training would stop if the validation loss
never improved again (`best epoch + patience`). The target moves out
whenever the validation loss improves, so the ETA is a running estimate,
not a promise.

## See also

[nsbi_parallel](https://neuralsbi.pedrodelima.com/reference/nsbi_parallel.md)
for running the simulator across cores.
