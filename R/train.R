#' Shared training engine for torch conditional density estimators
#'
#' All neural estimators (MDN, MAF, NSF) share one training loop so that
#' robustness features are implemented once: train/validation split, Adam,
#' minibatching, early stopping on validation loss, learning-rate decay on
#' plateau, gradient clipping, and best-of-`n_restarts` reinitialization.
#' The defaults (batch 200, lr 5e-4, 10% validation, patience 20, clip norm 5)
#' match Python `sbi`, so results are comparable across the two packages.
#'
#' @param build_net A zero-argument function returning a *fresh* torch module.
#'   Called once per restart so each restart gets new initial weights.
#' @param log_prob_fn `function(net, theta, x)` returning a length-`b` tensor of
#'   log densities for a `(b, p)` theta tensor and `(b, q)` x tensor.
#' @param theta,x Standardized training matrices.
#' @param n_restarts Train this many independently initialized networks and
#'   keep the one with the best validation loss.
#' @param clip_grad_norm Maximum gradient norm (set `Inf` to disable).
#' @param lr_patience,lr_factor,min_lr Reduce the learning rate by `lr_factor`
#'   after `lr_patience` epochs without validation improvement, down to
#'   `min_lr`.
#' @param device Torch device keyword: `"cpu"` (the default), `"cuda"`,
#'   `"mps"`, or `"gpu"`/`"auto"` (see `resolve_device()`, which is what turns
#'   this into an actual, available device). Training and validation tensors
#'   are created there, and the net is moved there right after `build_net()`,
#'   so the two never disagree the way they do under a bare
#'   `torch::with_device()`.
#' @return `list(net, best_val_loss, history, device)`, where `history` is a
#'   data frame of per-epoch train/validation losses for the winning restart
#'   and `device` is the resolved device (`"cpu"`, `"cuda"` or `"mps"`)
#'   training actually ran on.
#' @keywords internal
train_conditional_de <- function(build_net, log_prob_fn, theta, x,
                                 max_epochs = 2000L, batch_size = 200L,
                                 lr = 5e-4, validation_fraction = 0.1,
                                 patience = 20L, n_restarts = 1L,
                                 clip_grad_norm = 5,
                                 lr_patience = 10L, lr_factor = 0.5,
                                 min_lr = 1e-6,
                                 seed = NULL, verbose = FALSE,
                                 device = "cpu") {
  check_train_controls(max_epochs, batch_size, lr, validation_fraction,
                       patience, n_restarts, clip_grad_norm,
                       n = nrow(as_theta_matrix(theta)))
  # The bar spans every restart, so progress reporting is set up out here and
  # the loop itself lives in train_restarts().
  with_nsbi_progress(train_restarts(
    build_net, log_prob_fn, theta, x, max_epochs, batch_size, lr,
    validation_fraction, patience, n_restarts, clip_grad_norm, lr_patience,
    lr_factor, min_lr, seed, verbose, device
  ))
}

#' Validate the training controls
#'
#' Every one of these decides how the split, the batches or the restart loop
#' is built, and an unchecked value is reported by whichever base function
#' hits it first: `batch_size = 0` comes back as "invalid '(to - from)/by'",
#' `validation_fraction = 1` as "wrong sign in 'by' argument", and
#' `n_restarts = 0` as "Training failed: no restart produced a finite
#' validation loss", which blames training for an argument. [npe()] and [nle()]
#' call this before they simulate, so a typo does not cost the budget first.
#'
#' `n` is optional because that call happens before there are any rows. When it
#' is known, `validation_fraction` is checked against it: the requirement is
#' that both sides of the split come out non-empty, which the fraction alone
#' cannot decide.
#'
#' @inheritParams npe
#' @param n Number of training rows, or `NULL` when they do not exist yet.
#' @keywords internal
check_train_controls <- function(max_epochs, batch_size, lr,
                                 validation_fraction, patience, n_restarts,
                                 clip_grad_norm, n = NULL) {
  check_count(max_epochs, "max_epochs")
  check_count(batch_size, "batch_size")
  check_positive(lr, "lr")
  validation_fraction <- check_prob(validation_fraction, "validation_fraction")
  check_count(patience, "patience")
  check_count(n_restarts, "n_restarts", why = "since the best of them is kept")
  check_positive(clip_grad_norm, "clip_grad_norm", allow_inf = TRUE)

  if (!is.null(n)) {
    n_val <- max(1L, floor(validation_fraction * n))
    if (n - n_val < 1L) {
      need <- max(2L, ceiling(1 / (1 - validation_fraction)))
      stop(sprintf(paste0("`validation_fraction` of %s holds out %s of %d, ",
                          "leaving nothing to train on. At this fraction the ",
                          "estimator needs at least %s."),
                   format(validation_fraction), n_things(n_val, "row"), n,
                   n_things(as.integer(need), "row")),
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Restart loop behind `train_conditional_de()`
#'
#' Split out so `train_conditional_de()` can wrap the whole loop -- all
#' restarts, one progress bar -- in a progress-reporting context.
#'
#' `device` here is still the raw, unresolved keyword `train_conditional_de()`
#' was given (`"cpu"`, `"cuda"`, `"mps"`, `"gpu"` or `"auto"`); resolving it to
#' an actual, available device needs `torch` loaded (see `resolve_device()`),
#' and this is the first point that is guaranteed true -- `require_torch()` is
#' the line above. Doing it here rather than earlier in
#' `train_conditional_de()` keeps [check_train_controls()] (which needs no
#' torch at all) running first, so a bad `batch_size` is still reported before
#' an unavailable device is.
#' @keywords internal
train_restarts <- function(build_net, log_prob_fn, theta, x, max_epochs,
                           batch_size, lr, validation_fraction, patience,
                           n_restarts, clip_grad_norm, lr_patience, lr_factor,
                           min_lr, seed, verbose, device = "cpu") {
  require_torch()
  device <- resolve_device(device)
  if (!is.null(seed)) torch::torch_manual_seed(seed)
  theta <- as_theta_matrix(theta)
  x <- as_theta_matrix(x)
  n <- nrow(theta)

  # One split shared across restarts so validation losses are comparable.
  n_val <- max(1L, floor(validation_fraction * n))
  perm <- sample.int(n)
  val_idx <- perm[seq_len(n_val)]
  tr_idx <- perm[-seq_len(n_val)]

  # Built directly on `device`: torch_tensor() from an R matrix ignores
  # torch::with_device()'s default and always lands on CPU (the bug this
  # argument exists to work around), so the device has to be named here.
  tt <- torch::torch_tensor(theta, dtype = torch::torch_float(), device = device)
  xt <- torch::torch_tensor(x, dtype = torch::torch_float(), device = device)
  theta_tr <- tt[tr_idx, , drop = FALSE]; x_tr <- xt[tr_idx, , drop = FALSE]
  theta_val <- tt[val_idx, , drop = FALSE]; x_val <- xt[val_idx, , drop = FALSE]
  n_tr <- length(tr_idx)

  best <- list(net = NULL, val = Inf, history = NULL)

  # One progress step per epoch. The total is unknown up front (early stopping
  # decides it), so the bar targets the epoch training would stop at if the
  # validation loss never improved again, and revises that target upward every
  # time it does improve. See the "Training progress" section of ?nsbi_progress.
  max_epochs <- as.integer(max_epochs)
  n_restarts <- as.integer(n_restarts)
  epochs_done <- integer(0)   # epochs used by finished restarts
  p <- nsbi_progressor(steps = max_epochs * n_restarts, label = "Training")
  on.exit(p(0, done = TRUE), add = TRUE)

  for (restart in seq_len(n_restarts)) {
    net <- build_net()
    net$to(device = device)
    opt <- torch::optim_adam(net$parameters, lr = lr)
    scheduler <- torch::lr_reduce_on_plateau(opt, factor = lr_factor,
                                             patience = lr_patience,
                                             min_lr = min_lr)

    best_val <- Inf
    best_state <- NULL
    best_epoch <- 0L
    epochs_no_improve <- 0L
    hist_train <- numeric(0)
    hist_val <- numeric(0)

    for (epoch in seq_len(max_epochs)) {
      net$train()
      order <- sample.int(n_tr)
      starts <- seq(1L, n_tr, by = batch_size)
      epoch_loss <- 0
      for (s in starts) {
        idx <- order[s:min(s + batch_size - 1L, n_tr)]
        opt$zero_grad()
        lp <- log_prob_fn(net, theta_tr[idx, , drop = FALSE],
                          x_tr[idx, , drop = FALSE])
        loss <- -lp$mean()
        loss$backward()
        if (is.finite(clip_grad_norm)) {
          torch::nn_utils_clip_grad_norm_(net$parameters,
                                          max_norm = clip_grad_norm)
        }
        opt$step()
        epoch_loss <- epoch_loss + loss$item() * length(idx)
      }
      net$eval()
      val_loss <- torch::with_no_grad({
        (-log_prob_fn(net, theta_val, x_val)$mean())$item()
      })
      hist_train <- c(hist_train, epoch_loss / n_tr)
      hist_val <- c(hist_val, val_loss)

      if (is.finite(val_loss) && val_loss < best_val - 1e-4) {
        best_val <- val_loss
        best_state <- lapply(net$state_dict(), function(t) t$clone())
        best_epoch <- epoch
        epochs_no_improve <- 0L
      } else {
        epochs_no_improve <- epochs_no_improve + 1L
      }
      scheduler$step(val_loss)  # decay lr on validation plateau
      p(1, total = train_progress_total(epochs_done, best_epoch, patience,
                                        max_epochs, restart, n_restarts))
      if (verbose && (epoch %% 10L == 0L || epoch == 1L)) {
        cat(sprintf(
          "[train] restart %d epoch %d  val_loss=%.4f  best=%.4f\n",
          restart, epoch, val_loss, best_val))
      }
      if (epochs_no_improve >= patience) {
        verbose_cat(verbose, sprintf("[train] restart %d early stop at epoch %d\n",
                                     restart, epoch))
        break
      }
    }
    epochs_done <- c(epochs_done, epoch)
    if (!is.null(best_state)) net$load_state_dict(best_state)
    net$eval()

    if (best_val < best$val) {
      best <- list(
        net = net, val = best_val,
        history = data.frame(epoch = seq_along(hist_val),
                             train_loss = hist_train, val_loss = hist_val)
      )
    }
  }
  p(0, done = TRUE)

  if (is.null(best$net)) {
    stop("Training failed: no restart produced a finite validation loss.",
         call. = FALSE)
  }
  list(net = best$net, best_val_loss = best$val, history = best$history,
       device = device)
}

#' Projected epoch count for the training progress bar
#'
#' The run stops when `patience` epochs pass without a better validation loss,
#' so if the current best never improves it ends at `best_epoch + patience`.
#' That is the projection; restarts not yet started are budgeted at the mean
#' length of the ones already finished.
#' @keywords internal
train_progress_total <- function(epochs_done, best_epoch, patience, max_epochs,
                                 restart, n_restarts) {
  this <- min(max_epochs, best_epoch + patience)
  est <- if (length(epochs_done)) mean(epochs_done) else this
  sum(epochs_done) + this + (n_restarts - restart) * est
}
