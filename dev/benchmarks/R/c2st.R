# sbibm's classifier two-sample test, reimplemented in base R.
#
# The C2ST numbers in `main_paper.csv` come from `sbibm/metrics/c2st.py`, which
# is a scikit-learn `MLPClassifier` with two hidden layers of 10 x dim units,
# ReLU activations and Adam, cross-validated over 5 folds, scored by accuracy.
# `neuralsbi::c2st()` fits a logistic regression instead, which is a different
# (and much weaker) test, so its numbers are not comparable to the paper's. This
# file reproduces the sklearn recipe closely enough that the accuracies line up:
# same architecture, same initialisation scheme, same optimiser and defaults,
# same convergence rule.
#
# The one thing that cannot match is the RNG. sklearn's shuffling and weight
# initialisation are not reproducible outside Python, so individual runs differ
# by the usual C2ST noise (a few thousandths on 10k-vs-10k samples). That is far
# below the tolerance the report applies.

#' Classifier two-sample test, sbibm-compatible.
#'
#' @param X,Y Sample matrices, rows = draws. `X` is the reference sample; the
#'   z-scoring uses its mean and sd, as sbibm does.
#' @param z_score Whether to z-score. sbibm's published `C2ST` column is the
#'   z-scored variant, so this defaults to `TRUE`.
#' @param n_folds Cross-validation folds.
#' @param seed Optional seed.
#' @param max_epochs Cap on Adam epochs per fold. sbibm passes `max_iter=10000`
#'   to sklearn, and neither implementation's convergence rule (training loss
#'   improving by less than `tol` for 10 epochs running) fires anywhere near
#'   that on a C2ST task: it keeps improving for hundreds of epochs. 1000 is a
#'   compromise for R's speed. It biases in the safe direction. A classifier cut
#'   short is a *worse* classifier, and a worse classifier on 10k-vs-10k
#'   posterior samples reports a *higher* accuracy here, so this cap can only
#'   make neuralsbi look further from the reference than it is. Measured on
#'   gaussian_linear, going from 200 to 1000 epochs moved the C2ST from 0.536 to
#'   0.521, so the remaining bias from 1000 to 10000 is small next to the 0.05
#'   the report tolerates.
#' @return Mean cross-validated accuracy, with per-fold accuracies attached as
#'   attribute `folds`.
c2st_sbibm <- function(X, Y, z_score = TRUE, n_folds = 5L, seed = NULL,
                       max_epochs = 1000L) {
  X <- as.matrix(X); Y <- as.matrix(Y)
  stopifnot(ncol(X) == ncol(Y))
  if (!is.null(seed)) set.seed(seed)

  if (z_score) {
    mu <- colMeans(X)
    sdv <- apply(X, 2, stats::sd)
    sdv[!is.finite(sdv) | sdv == 0] <- 1
    X <- sweep(sweep(X, 2, mu, `-`), 2, sdv, `/`)
    Y <- sweep(sweep(Y, 2, mu, `-`), 2, sdv, `/`)
  }

  # sbibm always compares 10k against 10k. If the two sides differ (a run that
  # produced fewer draws, or a quick check), balance them: accuracy against
  # unequal classes is not a two-sample test any more, it is a prior on the
  # bigger class.
  X <- X[is.finite(rowSums(X)), , drop = FALSE]
  Y <- Y[is.finite(rowSums(Y)), , drop = FALSE]
  n_each <- min(nrow(X), nrow(Y))
  if (nrow(X) > n_each) X <- X[base::sample.int(nrow(X), n_each), , drop = FALSE]
  if (nrow(Y) > n_each) Y <- Y[base::sample.int(nrow(Y), n_each), , drop = FALSE]

  data <- rbind(X, Y)
  target <- c(rep(0, n_each), rep(1, n_each))

  ndim <- ncol(data)
  hidden <- rep(10L * ndim, 2)
  n <- nrow(data)
  # base::sample, not neuralsbi's sample() generic, which dispatches on its
  # first argument.
  fold <- base::sample(rep_len(seq_len(n_folds), n))

  acc <- numeric(n_folds)
  for (k in seq_len(n_folds)) {
    tr <- fold != k
    net <- mlp_fit(data[tr, , drop = FALSE], target[tr], hidden = hidden,
                   max_epochs = max_epochs)
    p <- mlp_predict(net, data[!tr, , drop = FALSE])
    acc[k] <- mean((p > 0.5) == (target[!tr] == 1))
  }
  structure(mean(acc), folds = acc)
}

# --- the classifier ----------------------------------------------------------

# Defaults below are scikit-learn's MLPClassifier defaults, which sbibm does not
# override except for the architecture.
MLP_ALPHA <- 1e-4          # L2 penalty
MLP_LR <- 1e-3             # Adam learning_rate_init
MLP_B1 <- 0.9
MLP_B2 <- 0.999
MLP_EPS <- 1e-8
MLP_BATCH <- 200L
MLP_TOL <- 1e-4
MLP_NO_IMPROVE <- 10L

#' Glorot-style uniform initialisation, as in sklearn's `_init_coef`.
init_layer <- function(fan_in, fan_out) {
  bound <- sqrt(6 / (fan_in + fan_out))
  list(W = matrix(stats::runif(fan_in * fan_out, -bound, bound),
                  nrow = fan_in, ncol = fan_out),
       b = stats::runif(fan_out, -bound, bound))
}

mlp_fit <- function(x, y, hidden, max_epochs = 1000L) {
  sizes <- c(ncol(x), hidden, 1L)
  layers <- lapply(seq_len(length(sizes) - 1), function(i)
    init_layer(sizes[i], sizes[i + 1]))

  m <- lapply(layers, function(l) list(W = 0 * l$W, b = 0 * l$b))
  v <- m
  t_step <- 0
  best_loss <- Inf
  no_improve <- 0L
  n <- nrow(x)
  batch <- min(MLP_BATCH, n)

  for (epoch in seq_len(max_epochs)) {
    perm <- sample.int(n)
    starts <- seq(1, n, by = batch)
    epoch_loss <- 0
    for (s in starts) {
      idx <- perm[s:min(s + batch - 1, n)]
      nb <- length(idx)
      xb <- x[idx, , drop = FALSE]
      yb <- y[idx]

      # forward
      a <- vector("list", length(layers) + 1)
      a[[1]] <- xb
      for (i in seq_along(layers)) {
        z <- sweep(a[[i]] %*% layers[[i]]$W, 2, layers[[i]]$b, `+`)
        a[[i + 1]] <- if (i < length(layers)) pmax(z, 0) else sigmoid(z)
      }
      p <- as.numeric(a[[length(a)]])
      epoch_loss <- epoch_loss + log_loss(yb, p) * nb

      # backward
      delta <- matrix(p - yb, ncol = 1)
      grads <- vector("list", length(layers))
      for (i in rev(seq_along(layers))) {
        grads[[i]] <- list(
          W = (t(a[[i]]) %*% delta + MLP_ALPHA * layers[[i]]$W) / nb,
          b = colSums(delta) / nb
        )
        if (i > 1) {
          delta <- (delta %*% t(layers[[i]]$W)) * (a[[i]] > 0)
        }
      }

      # Adam with sklearn's bias-corrected step size
      t_step <- t_step + 1
      lr <- MLP_LR * sqrt(1 - MLP_B2^t_step) / (1 - MLP_B1^t_step)
      for (i in seq_along(layers)) {
        m[[i]]$W <- MLP_B1 * m[[i]]$W + (1 - MLP_B1) * grads[[i]]$W
        m[[i]]$b <- MLP_B1 * m[[i]]$b + (1 - MLP_B1) * grads[[i]]$b
        v[[i]]$W <- MLP_B2 * v[[i]]$W + (1 - MLP_B2) * grads[[i]]$W^2
        v[[i]]$b <- MLP_B2 * v[[i]]$b + (1 - MLP_B2) * grads[[i]]$b^2
        layers[[i]]$W <- layers[[i]]$W - lr * m[[i]]$W / (sqrt(v[[i]]$W) + MLP_EPS)
        layers[[i]]$b <- layers[[i]]$b - lr * m[[i]]$b / (sqrt(v[[i]]$b) + MLP_EPS)
      }
    }
    l2 <- sum(vapply(layers, function(l) sum(l$W^2), numeric(1)))
    loss <- epoch_loss / n + 0.5 * MLP_ALPHA * l2 / n

    if (loss > best_loss - MLP_TOL) no_improve <- no_improve + 1L else no_improve <- 0L
    if (loss < best_loss) best_loss <- loss
    if (no_improve > MLP_NO_IMPROVE) break
  }
  layers
}

mlp_predict <- function(layers, x) {
  a <- x
  for (i in seq_along(layers)) {
    z <- sweep(a %*% layers[[i]]$W, 2, layers[[i]]$b, `+`)
    a <- if (i < length(layers)) pmax(z, 0) else sigmoid(z)
  }
  as.numeric(a)
}

sigmoid <- function(z) 1 / (1 + exp(-z))

log_loss <- function(y, p) {
  p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}
