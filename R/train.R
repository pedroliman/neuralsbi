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
#' @return `list(net, best_val_loss, history)` where `history` is a data frame
#'   of per-epoch train/validation losses for the winning restart.
#' @keywords internal
train_conditional_de <- function(build_net, log_prob_fn, theta, x,
                                 max_epochs = 2000L, batch_size = 200L,
                                 lr = 5e-4, validation_fraction = 0.1,
                                 patience = 20L, n_restarts = 1L,
                                 clip_grad_norm = 5,
                                 lr_patience = 10L, lr_factor = 0.5,
                                 min_lr = 1e-6,
                                 seed = NULL, verbose = FALSE) {
  # The bar spans every restart, so progress reporting is set up out here and
  # the loop itself lives in train_restarts().
  with_nsbi_progress(train_restarts(
    build_net, log_prob_fn, theta, x, max_epochs, batch_size, lr,
    validation_fraction, patience, n_restarts, clip_grad_norm, lr_patience,
    lr_factor, min_lr, seed, verbose
  ))
}

#' Restart loop behind `train_conditional_de()`
#'
#' Split out so `train_conditional_de()` can wrap the whole loop -- all
#' restarts, one progress bar -- in a progress-reporting context.
#' @keywords internal
train_restarts <- function(build_net, log_prob_fn, theta, x, max_epochs,
                           batch_size, lr, validation_fraction, patience,
                           n_restarts, clip_grad_norm, lr_patience, lr_factor,
                           min_lr, seed, verbose) {
  require_torch()
  if (!is.null(seed)) torch::torch_manual_seed(seed)
  theta <- as_theta_matrix(theta)
  x <- as_theta_matrix(x)
  n <- nrow(theta)

  # One split shared across restarts so validation losses are comparable.
  n_val <- max(1L, floor(validation_fraction * n))
  perm <- sample.int(n)
  val_idx <- perm[seq_len(n_val)]
  tr_idx <- perm[-seq_len(n_val)]

  tt <- torch::torch_tensor(theta, dtype = torch::torch_float())
  xt <- torch::torch_tensor(x, dtype = torch::torch_float())
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
        verbose_cat(TRUE, sprintf(
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
  list(net = best$net, best_val_loss = best$val, history = best$history)
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
