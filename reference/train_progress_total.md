# Projected epoch count for the training progress bar

The run stops when `patience` epochs pass without a better validation
loss, so if the current best never improves it ends at
`best_epoch + patience`. That is the projection; restarts not yet
started are budgeted at the mean length of the ones already finished.

## Usage

``` r
train_progress_total(
  epochs_done,
  best_epoch,
  patience,
  max_epochs,
  restart,
  n_restarts
)
```
