#' Classifier two-sample test (C2ST)
#'
#' Trains a classifier to tell the draws in `x` apart from the draws in `y` and
#' reports its cross-validated test accuracy. An accuracy near 0.5 means the two
#' sample sets are indistinguishable (good); near 1.0 means they differ. This is
#' the headline metric of the `sbibm` benchmark suite, used there to score an
#' approximate posterior against reference draws.
#'
#' The defaults reproduce the procedure in `sbibm/metrics/c2st.py`, so a number
#' from here is comparable with a published one: z-score both sample sets using
#' the mean and standard deviation of `x`, fit a two-hidden-layer ReLU network
#' of `10 * d` units per layer (`d` = number of columns) with Adam, L2 penalty
#' `1e-4`, minibatches of 200 and a learning rate of `1e-3`, and average the
#' accuracy over 5 shuffled cross-validation folds. Training stops when the
#' epoch loss fails to improve by `1e-4` for 10 epochs in a row, which is what
#' ends the fit long before `max_epochs`. `sbibm` passes reference draws first,
#' so pass the reference as `x`: that is the set whose moments set the scale.
#'
#' The network is written in base R rather than delegating to torch, because
#' C2ST is a diagnostic and must work in a session with no torch installed.
#' Budget about six seconds for the 10000-against-10000 comparison `sbibm` runs
#' in two dimensions, and more as `d` grows, since the hidden layers grow with
#' it. Where that is too slow, lower `hidden` or `max_epochs`, or take the
#' logistic screen below.
#'
#' `classifier = "logistic"` swaps the network for cross-validated logistic
#' regression. It is fast and deterministic, and it is what this function used
#' before it was aligned with `sbibm`, but it is linear: it sees a shift in
#' location and is close to blind to two sample sets that share a mean and
#' differ in spread or in the shape of their dependence. Use it as a cheap
#' screen, and report the MLP number.
#'
#' Unequal sample sizes are balanced by subsampling the larger set, which
#' `sbibm` does not do because it always compares 10000 draws against 10000.
#' Accuracy against unbalanced classes is not a two-sample test: 8000 draws
#' against 2000 identical ones scores 0.8 for a classifier that has learned
#' nothing except to always answer with the bigger class.
#'
#' @param x,y Matrices of samples (rows = draws, cols = dimensions). The two
#'   must be the same width, since they are draws of the same quantity. Row
#'   counts need not match; the larger set is subsampled down to the smaller.
#'   Pass reference draws as `x`, following `sbibm`.
#' @param n_folds Number of cross-validation folds. At least 2, and fewer than
#'   the number of draws in the smaller sample set.
#' @param seed Optional seed. Fixes the fold split, the subsampling, the network
#'   initialization and the minibatch order.
#' @param classifier `"mlp"` for the `sbibm` network, `"logistic"` for
#'   logistic regression.
#' @param z_score Standardize both sample sets by the mean and standard
#'   deviation of `x` before training. On by default, as in `sbibm`.
#' @param noise_scale Standard deviation of Gaussian noise added to both sample
#'   sets after z-scoring. `NULL`, the default, adds none. Set it when the draws
#'   are discrete or lie on a lower-dimensional set, where a classifier can
#'   separate the two sides on an artefact of representation.
#' @param hidden Widths of the hidden layers, one entry per layer. `NULL`, the
#'   default, is `sbibm`'s two layers of `10 * d` units.
#' @param max_epochs Cap on training epochs per fold. A guard, not a budget:
#'   the no-improvement rule normally stops training first.
#' @return A list with the mean cross-validated accuracy and ROC AUC, the
#'   per-fold values of each, the number of draws per class, and a one-line
#'   reading of the accuracy.
#' @references Lopez-Paz, D. and Oquab, M. (2017). Revisiting classifier
#'   two-sample tests. *ICLR*. \doi{10.48550/arXiv.1610.06545}
#'
#'   Lueckmann, J.-M., Boelts, J., Greenberg, D. S., Goncalves, P. J. and
#'   Macke, J. H. (2021). Benchmarking simulation-based inference. *AISTATS*.
#'   \doi{10.48550/arXiv.2101.04653}
#' @examples
#' a <- matrix(rnorm(400), ncol = 2)
#' b <- matrix(rnorm(400), ncol = 2)
#' c2st(a, b, classifier = "logistic", seed = 1)$accuracy
#' @export
c2st <- function(x, y, n_folds = 5L, seed = NULL,
                 classifier = c("mlp", "logistic"), z_score = TRUE,
                 noise_scale = NULL, hidden = NULL, max_epochs = 10000L) {
  classifier <- match.arg(classifier)
  n_folds <- check_count(n_folds, "n_folds", min = 2L,
                         why = "since each fold is scored by a fit on the rest")
  max_epochs <- check_count(max_epochs, "max_epochs")
  if (!is.null(noise_scale)) noise_scale <- check_positive(noise_scale, "noise_scale")
  if (!is.null(seed)) set.seed(seed)
  # A row is one draw here, so a bare vector is a column of 1-D draws rather
  # than check_matrix()'s single row. That is the pre-computed (theta, x) rule
  # in npe(), and it is why the type check and the reshape are separate calls.
  x <- as_theta_matrix(check_numeric(x, "x"))
  y <- as_theta_matrix(check_numeric(y, "y"))
  if (ncol(x) != ncol(y)) {
    # rbind() is where this lands otherwise, and it complains about the number
    # of columns of a matrix it was handed rather than about x or y.
    stop(sprintf(paste0("`x` has %s and `y` has %s. c2st() compares two sets ",
                        "of draws of the same quantity, so both must be the ",
                        "same width (rows = draws, columns = dimensions)."),
                 n_things(ncol(x), "column"), n_things(ncol(y), "column")),
         call. = FALSE)
  }
  n_each <- min(nrow(x), nrow(y))
  if (n_folds >= n_each) {
    # Folds are cut over both sample sets together, so n_folds up to n_each
    # still leaves each training fit some draws from each class. Past that the
    # test folds thin out and empty ones score NaN, which propagates to the
    # accuracy with nothing said about why.
    stop(sprintf(paste0("`n_folds` is %d, but the smaller sample set has only ",
                        "%s. It must be fewer, so that every fold has draws to ",
                        "test on."),
                 n_folds, n_things(n_each, "draw")),
         call. = FALSE)
  }
  d <- ncol(x)
  hidden <- if (is.null(hidden)) rep(10L * d, 2L) else
    check_counts(hidden, "hidden", what = "hidden units")
  # base::sample.int, not the package's sample() generic, which dispatches on
  # its first argument.
  if (nrow(x) > n_each) x <- x[base::sample.int(nrow(x), n_each), , drop = FALSE]
  if (nrow(y) > n_each) y <- y[base::sample.int(nrow(y), n_each), , drop = FALSE]
  if (isTRUE(z_score)) {
    # sbibm scales both sides by the moments of X alone, not of the pool. The
    # reference set fixes the units, so the number does not move when the
    # sample set being scored is the one that is wrong.
    std <- fit_standardizer(x)
    x <- apply_standardizer(std, x)
    y <- apply_standardizer(std, y)
  }
  if (!is.null(noise_scale)) {
    x <- x + matrix(stats::rnorm(length(x), sd = noise_scale), nrow(x))
    y <- y + matrix(stats::rnorm(length(y), sd = noise_scale), nrow(y))
  }
  data <- rbind(x, y)
  # sbibm labels X 0 and Y 1. Accuracy does not care which way round, but the
  # AUC does, so keep the convention.
  label <- c(rep(0, nrow(x)), rep(1, nrow(y)))
  n <- nrow(data)
  # base::sample again: sample() is an S3 generic in this package.
  fold <- base::sample(rep_len(seq_len(n_folds), n))
  accs <- numeric(n_folds)
  aucs <- numeric(n_folds)
  for (k in seq_len(n_folds)) {
    tr <- fold != k
    te <- !tr
    prob <- switch(
      classifier,
      mlp = c2st_mlp_prob(data[tr, , drop = FALSE], label[tr],
                          data[te, , drop = FALSE], hidden = hidden,
                          max_epochs = max_epochs),
      logistic = c2st_logistic_prob(data[tr, , drop = FALSE], label[tr],
                                    data[te, , drop = FALSE])
    )
    accs[k] <- mean((prob > 0.5) == (label[te] == 1))
    aucs[k] <- roc_auc(prob, label[te])
  }
  list(accuracy = mean(accs), auc = mean(aucs),
       fold_accuracy = accs, fold_auc = aucs,
       n = n_each, classifier = classifier,
       interpretation = if (mean(accs) < 0.55)
         "indistinguishable (good)" else "distinguishable")
}

#' Cross-validated logistic regression for [c2st()]
#'
#' @param x_train,y_train Training draws and their 0/1 labels.
#' @param x_test Draws to score.
#' @return Predicted probability of class 1 for each row of `x_test`.
#' @keywords internal
c2st_logistic_prob <- function(x_train, y_train, x_test) {
  df <- data.frame(y = y_train, x_train)
  fit <- suppressWarnings(stats::glm(y ~ ., data = df, family = stats::binomial()))
  stats::predict(fit, newdata = data.frame(x_test), type = "response")
}

#' `sbibm`'s C2ST classifier
#'
#' A two-hidden-layer ReLU network trained by Adam on the binary log loss,
#' matching the `MLPClassifier` settings `sbibm` uses: Glorot-uniform
#' initialization, L2 penalty `alpha`, minibatches of 200 drawn in a fresh
#' random order each epoch, Adam with bias-corrected step size, and a stop once
#' the epoch loss has failed to improve on the best by more than `tol` for
#' `n_iter_no_change` epochs in a row.
#'
#' Written out in base R on purpose. C2ST is a diagnostic, and the diagnostics
#' have to run in a session where torch was never installed.
#'
#' @param x_train,y_train Training draws and their 0/1 labels.
#' @param x_test Draws to score.
#' @param hidden Widths of the hidden layers, one entry per layer.
#' @param alpha L2 penalty, scaled by batch size the way `scikit-learn` scales
#'   it.
#' @param lr Adam step size.
#' @param batch_size Minibatch size, capped at the number of training draws.
#' @param max_epochs Cap on epochs.
#' @param tol Smallest loss improvement that counts as progress.
#' @param n_iter_no_change Epochs without progress before training stops.
#' @return Predicted probability of class 1 for each row of `x_test`.
#' @keywords internal
c2st_mlp_prob <- function(x_train, y_train, x_test, hidden,
                          alpha = 1e-4, lr = 1e-3, batch_size = 200L,
                          max_epochs = 10000L, tol = 1e-4,
                          n_iter_no_change = 10L) {
  net <- c2st_mlp_train(x_train, y_train, hidden = hidden, alpha = alpha,
                        lr = lr, batch_size = batch_size,
                        max_epochs = max_epochs, tol = tol,
                        n_iter_no_change = n_iter_no_change)
  c2st_mlp_forward(net, x_test)
}

#' Add a bias row-wise
#'
#' `sweep()` and `t(t(z) + b)` both copy more than they need to; this is the
#' same thing written as one recycled vector, and the inner loop of the
#' network runs it a few hundred thousand times.
#'
#' @param z An `n x p` matrix.
#' @param b A length-`p` bias vector.
#' @keywords internal
c2st_add_bias <- function(z, b) z + rep(b, each = nrow(z))

#' Forward pass of the [c2st()] network
#'
#' @param net Weights and biases from [c2st_mlp_train()].
#' @param x Draws to score.
#' @return Predicted probability of class 1 for each row of `x`.
#' @keywords internal
c2st_mlp_forward <- function(net, x) {
  n_layers <- length(net$w)
  a <- x
  for (i in seq_len(n_layers)) {
    z <- c2st_add_bias(a %*% net$w[[i]], net$b[[i]])
    a <- if (i < n_layers) pmax(z, 0) else stats::plogis(z)
  }
  a[, 1L]
}

#' Train the [c2st()] network
#'
#' @inheritParams c2st_mlp_prob
#' @return A list of weight matrices `w`, bias vectors `b`, and the number of
#'   epochs run.
#' @keywords internal
c2st_mlp_train <- function(x_train, y_train, hidden, alpha = 1e-4, lr = 1e-3,
                           batch_size = 200L, max_epochs = 10000L, tol = 1e-4,
                           n_iter_no_change = 10L) {
  sizes <- c(ncol(x_train), hidden, 1L)
  n_layers <- length(sizes) - 1L
  w <- vector("list", n_layers)
  b <- vector("list", n_layers)
  for (i in seq_len(n_layers)) {
    # Glorot uniform, as scikit-learn initializes a ReLU layer.
    bound <- sqrt(6 / (sizes[i] + sizes[i + 1L]))
    w[[i]] <- matrix(stats::runif(sizes[i] * sizes[i + 1L], -bound, bound),
                     sizes[i], sizes[i + 1L])
    b[[i]] <- stats::runif(sizes[i + 1L], -bound, bound)
  }
  zeros_like <- function(p) lapply(p, function(q) q * 0)
  mw <- zeros_like(w); vw <- zeros_like(w)
  mb <- zeros_like(b); vb <- zeros_like(b)
  beta1 <- 0.9; beta2 <- 0.999; adam_eps <- 1e-8
  n <- nrow(x_train)
  batch_size <- min(as.integer(batch_size), n)
  clip <- .Machine$double.eps
  step <- 0L
  best_loss <- Inf
  no_improvement <- 0L
  epochs <- 0L
  for (epoch in seq_len(max_epochs)) {
    epochs <- epoch
    ord <- base::sample.int(n)
    running <- 0
    start <- 1L
    while (start <= n) {
      idx <- ord[start:min(start + batch_size - 1L, n)]
      start <- start + batch_size
      nb <- length(idx)
      xb <- x_train[idx, , drop = FALSE]
      yb <- y_train[idx]

      a <- vector("list", n_layers + 1L)
      a[[1L]] <- xb
      for (i in seq_len(n_layers)) {
        z <- c2st_add_bias(a[[i]] %*% w[[i]], b[[i]])
        a[[i + 1L]] <- if (i < n_layers) pmax(z, 0) else stats::plogis(z)
      }
      out <- a[[n_layers + 1L]][, 1L]
      p <- pmin(pmax(out, clip), 1 - clip)
      loss <- -sum(yb * log(p) + (1 - yb) * log1p(-p)) / nb +
        0.5 * alpha * sum(vapply(w, function(m) sum(m * m), 0)) / nb
      running <- running + loss * nb

      delta <- matrix(out - yb, ncol = 1L)
      for (i in rev(seq_len(n_layers))) {
        gw <- (crossprod(a[[i]], delta) + alpha * w[[i]]) / nb
        gb <- colMeans(delta)
        if (i > 1L) delta <- tcrossprod(delta, w[[i]]) * (a[[i]] > 0)
        step_i <- step + 1L
        # Adam's bias correction folded into the step size, as scikit-learn
        # writes it, so the first steps are not shrunk toward zero.
        lr_t <- lr * sqrt(1 - beta2^step_i) / (1 - beta1^step_i)
        mw[[i]] <- beta1 * mw[[i]] + (1 - beta1) * gw
        vw[[i]] <- beta2 * vw[[i]] + (1 - beta2) * gw * gw
        w[[i]] <- w[[i]] - lr_t * mw[[i]] / (sqrt(vw[[i]]) + adam_eps)
        mb[[i]] <- beta1 * mb[[i]] + (1 - beta1) * gb
        vb[[i]] <- beta2 * vb[[i]] + (1 - beta2) * gb * gb
        b[[i]] <- b[[i]] - lr_t * mb[[i]] / (sqrt(vb[[i]]) + adam_eps)
      }
      step <- step + 1L
    }
    epoch_loss <- running / n
    no_improvement <- if (epoch_loss > best_loss - tol) no_improvement + 1L else 0L
    if (epoch_loss < best_loss) best_loss <- epoch_loss
    if (no_improvement > n_iter_no_change) break
  }
  list(w = w, b = b, epochs = epochs)
}

#' Area under the ROC curve
#'
#' The rank form of the Mann-Whitney statistic, which needs no threshold sweep
#' and handles ties the way `scikit-learn`'s `roc_auc_score` does.
#'
#' @param prob Predicted probability of class 1.
#' @param label 0/1 labels.
#' @keywords internal
roc_auc <- function(prob, label) {
  pos <- label == 1
  n1 <- sum(pos)
  n0 <- length(label) - n1
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(prob)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
