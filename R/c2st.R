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
#' accuracy over 5 shuffled cross-validation folds, stratified by class the way
#' `sklearn`'s default `StratifiedKFold` is (`sbibm` reaches it via
#' `cross_val_score()`). Training stops when the epoch loss fails to improve
#' by `1e-4` for 10 epochs in a row, which is what
#' ends the fit long before `max_epochs`. `sbibm` passes reference draws first,
#' so pass the reference as `x`: that is the set whose moments set the scale.
#'
#' The network trains on torch, like every other neural piece of this package,
#' so `classifier = "mlp"` needs torch installed and errors without it. That is
#' the only part of the diagnostics that does.
#'
#' Budget for the fit. `sbibm`'s own comparison, 10000 draws against 10000 in
#' two dimensions, takes about 25 seconds on one CPU core; at `d = 5` it is
#' several minutes, because the hidden layers are `10 * d` wide and the loss
#' takes longer to flatten. Five folds of an all but unbounded epoch budget is
#' what the `sbibm` recipe asks for, and that is what it costs. `hidden` and
#' `max_epochs` cut it down, and `device` moves the fit to a GPU.
#'
#' `classifier = "logistic"` swaps the network for cross-validated logistic
#' regression. It needs no torch, it answers in milliseconds, and it is what
#' this function used before it was aligned with `sbibm`. It is also linear: it
#' sees a shift in location and is close to blind to two sample sets that share
#' a mean and differ in spread or in the shape of their dependence. Use it as a
#' cheap screen, or in a session with no torch, and report the MLP number.
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
#' @param device Torch device to train the network on: `"cpu"` (the default),
#'   `"cuda"`, `"mps"`, or `"gpu"`/`"auto"` for whichever is available. Ignored
#'   by `classifier = "logistic"`.
#' @return A list with the mean cross-validated accuracy and ROC AUC, the
#'   per-fold values of each, the number of draws per class, the classifier
#'   used, and a one-line reading of the accuracy.
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
                 noise_scale = NULL, hidden = NULL, max_epochs = 10000L,
                 device = "cpu") {
  classifier <- match.arg(classifier)
  n_folds <- check_count(n_folds, "n_folds", min = 2L,
                         why = "since each fold is scored by a fit on the rest")
  max_epochs <- check_count(max_epochs, "max_epochs")
  check_device_arg(device)
  if (!is.null(hidden)) {
    hidden <- check_counts(hidden, "hidden", what = "hidden units")
  }
  if (!is.null(noise_scale)) noise_scale <- check_positive(noise_scale, "noise_scale")
  if (!is.null(seed)) set.seed(seed)
  # A row is one draw here, so a bare vector is a column of 1-D draws rather
  # than check_matrix()'s single row. That is the pre-computed (theta, x) rule
  # in npe(), and it is why the type check and the reshape are separate calls.
  # check_numeric() only enforces type, so an NA/NaN/Inf entry used to
  # standardize via fit_standardizer()/apply_standardizer() into a whole
  # column of NA (colMeans()/sd() carry no na.rm), corrupting every draw in
  # that column instead of just the offending row (#222). Checked explicitly
  # here, as surrogate_score() (#202) and mcmc_log_prob() (#208/#209) already
  # do for the same failure mode.
  x <- check_numeric(x, "x")
  check_finite(x, "x")
  y <- check_numeric(y, "y")
  check_finite(y, "y")
  x <- as_theta_matrix(x)
  y <- as_theta_matrix(y)
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
    # Folds are cut within x and within y separately (see
    # c2st_stratified_folds()), so n_folds up to n_each still leaves every
    # fold's test set at least one draw of each class. Past that the per-class
    # folds thin out and an empty one scores NA (roc_auc()'s n1 == 0 or
    # n0 == 0 guard), which propagates through mean(aucs) with nothing said
    # about why.
    stop(sprintf(paste0("`n_folds` is %d, but the smaller sample set has only ",
                        "%s. It must be fewer, so that every fold has draws to ",
                        "test on."),
                 n_folds, n_things(n_each, "draw")),
         call. = FALSE)
  }
  if (is.null(hidden)) hidden <- rep(10L * ncol(x), 2L)
  if (identical(classifier, "mlp")) {
    # Last, so that a bad argument or a mismatched pair of sample sets is
    # reported as itself rather than as a missing torch install.
    require_torch(
      what = "c2st()'s MLP classifier, the one sbibm uses,",
      alternative = paste("Alternatively use classifier = \"logistic\" for a",
                          "torch-free two-sample test. It is linear, so it",
                          "sees a shift in location but not a difference in",
                          "spread; see ?c2st.")
    )
    device <- resolve_device(device)
    # Two random streams to fix: R's decides the folds, the subsampling and the
    # batch order, torch's decides the network's starting weights.
    if (!is.null(seed)) {
      old_torch_rng <- set_torch_seed(seed)
      on.exit(torch::torch_set_rng_state(old_torch_rng), add = TRUE)
    }
  }
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
  fold <- c2st_stratified_folds(nrow(x), nrow(y), n_folds)
  accs <- numeric(n_folds)
  aucs <- numeric(n_folds)
  for (k in seq_len(n_folds)) {
    tr <- fold != k
    te <- !tr
    prob <- switch(
      classifier,
      mlp = c2st_mlp_prob(data[tr, , drop = FALSE], label[tr],
                          data[te, , drop = FALSE], hidden = hidden,
                          max_epochs = max_epochs, device = device),
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

#' Class-stratified cross-validation folds for [c2st()]
#'
#' Assigns fold labels to the x-rows and the y-rows of `data <- rbind(x, y)`
#' separately, then concatenates the two, rather than shuffling one fold
#' assignment over both classes at once. Unstratified shuffling can hand a
#' fold's test set zero rows of one class whenever a class is small or the
#' shuffle is unlucky, and [roc_auc()] returns `NA_real_` for that fold, so
#' `c2st()$auc` silently comes back `NA` (#269). This is what `sklearn`'s
#' `StratifiedKFold` guarantees for a two-class target, which is what `sbibm`
#' gets from `cross_val_score()`'s default cv for a classification estimator.
#'
#' `n_folds < min(n_x, n_y)` (enforced by the caller) keeps every per-class
#' `rep_len()` split at one row or more per fold, so every fold's test set
#' still has at least one draw of each class.
#'
#' @param n_x,n_y Number of draws in each class; x-rows come first in `data`.
#' @param n_folds Number of folds.
#' @return A length-`n_x + n_y` vector of fold indices in `seq_len(n_folds)`.
#' @keywords internal
c2st_stratified_folds <- function(n_x, n_y, n_folds) {
  # base::sample, not the package's sample() generic, which dispatches on its
  # first argument.
  c(base::sample(rep_len(seq_len(n_folds), n_x)),
    base::sample(rep_len(seq_len(n_folds), n_y)))
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

#' Build `sbibm`'s C2ST classifier as a torch module
#'
#' A plain `linear/relu` stack ending in one logit, the same trunk
#' [mlp_layers()] builds for every other estimator in the package. The one
#' thing worth doing by hand is the initialization: `scikit-learn` draws both
#' weights and biases from Glorot uniform, `torch` uses Kaiming uniform on the
#' weight and a fan-in bound on the bias, and the fit has to start where
#' `MLPClassifier`'s would for the numbers to line up.
#'
#' @param d Number of columns in the draws.
#' @param hidden Widths of the hidden layers, one entry per layer.
#' @keywords internal
c2st_mlp_module <- function(d, hidden) {
  torch::nn_module(
    classname = "nsbi_c2st_net",
    initialize = function() {
      layers <- mlp_layers(c(d, hidden))
      layers[[length(layers) + 1L]] <-
        torch::nn_linear(hidden[length(hidden)], 1L)
      for (layer in layers) {
        if (!inherits(layer, "nn_linear")) next
        bound <- sqrt(6 / (layer$in_features + layer$out_features))
        torch::nn_init_uniform_(layer$weight, -bound, bound)
        torch::nn_init_uniform_(layer$bias, -bound, bound)
      }
      self$trunk <- do.call(torch::nn_sequential, layers)
    },
    forward = function(x) self$trunk(x)$squeeze(2)
  )
}

#' `sbibm`'s C2ST classifier, trained on one fold
#'
#' Adam on the binary log loss under the `MLPClassifier` settings `sbibm` uses:
#' minibatches of 200 in a fresh random order each epoch, learning rate `1e-3`,
#' and a stop once the epoch loss has failed to improve on the best by more
#' than `tol` for `n_iter_no_change` epochs in a row.
#'
#' `scikit-learn`'s `alpha` penalizes the weights and leaves the biases alone,
#' adding `alpha * W / batch` to the gradient. That is what `weight_decay`
#' does, so the optimizer carries it in two parameter groups rather than the
#' loss carrying it as a graph node: the penalty is also part of the loss the
#' stopping rule reads, but evaluating it once an epoch instead of once a batch
#' costs nothing and saves half the running time of the loop.
#'
#' @param x_train,y_train Training draws and their 0/1 labels.
#' @param x_test Draws to score.
#' @param hidden Widths of the hidden layers, one entry per layer.
#' @param alpha L2 penalty on the weights.
#' @param lr Adam step size.
#' @param batch_size Minibatch size, capped at the number of training draws.
#' @param max_epochs Cap on epochs.
#' @param tol Smallest loss improvement that counts as progress.
#' @param n_iter_no_change Epochs without progress before training stops.
#' @param device Resolved torch device to train on.
#' @return Predicted probability of class 1 for each row of `x_test`.
#' @keywords internal
c2st_mlp_prob <- function(x_train, y_train, x_test, hidden, alpha = 1e-4,
                          lr = 1e-3, batch_size = 200L, max_epochs = 10000L,
                          tol = 1e-4, n_iter_no_change = 10L, device = "cpu") {
  net <- c2st_mlp_module(ncol(x_train), hidden)()
  net$to(device = device)
  n <- nrow(x_train)
  batch_size <- min(as.integer(batch_size), n)

  # A weight is the two-dimensional parameter of a linear layer, a bias the
  # one-dimensional one, and only the weights are penalized.
  is_weight <- function(p) length(dim(p)) == 2L
  weights <- Filter(is_weight, net$parameters)
  opt <- c2st_adam(
    list(list(params = weights, weight_decay = alpha / batch_size),
         list(params = Filter(Negate(is_weight), net$parameters),
              weight_decay = 0)),
    lr = lr)
  # Named device, not with_device(): torch_tensor() from an R matrix always
  # lands on the CPU otherwise, the same wrinkle train_restarts() works around.
  xt <- torch::torch_tensor(x_train, dtype = torch::torch_float(),
                            device = device)
  yt <- torch::torch_tensor(y_train, dtype = torch::torch_float(),
                            device = device)

  best_loss <- Inf
  no_improvement <- 0L
  for (epoch in seq_len(max_epochs)) {
    net$train()
    running <- 0
    # sample.int() and minibatches(), so R's seed governs the batch order the
    # way it governs every other training loop in the package.
    for (idx in minibatches(base::sample.int(n), batch_size)) {
      opt$zero_grad()
      logit <- net(xt[idx, , drop = FALSE])
      loss <- torch::nnf_binary_cross_entropy_with_logits(logit, yt[idx])
      loss$backward()
      opt$step()
      running <- running + loss$item() * length(idx)
    }
    # The penalty term scikit-learn folds into the loss it watches, on the
    # weights as they end the epoch.
    penalty <- torch::with_no_grad(
      0.5 * alpha * sum(vapply(weights, function(w) {
        as.numeric(torch::torch_sum(w * w)$item())
      }, 0)) / batch_size)
    epoch_loss <- running / n + penalty
    no_improvement <- if (epoch_loss > best_loss - tol) no_improvement + 1L else 0L
    if (epoch_loss < best_loss) best_loss <- epoch_loss
    if (no_improvement > n_iter_no_change) break
  }

  net$eval()
  torch::with_no_grad({
    xte <- torch::torch_tensor(x_test, dtype = torch::torch_float(),
                               device = device)
    as.numeric(torch::torch_sigmoid(net(xte))$cpu())
  })
}

#' Adam for [c2st_mlp_prob()], preferring torch's C++ implementation
#'
#' `optim_ignite_adam()` takes the whole step in C++ where `optim_adam()` loops
#' over the parameters in R. The two produce identical trajectories, and this
#' loop is thousands of steps over six small tensors, so the R-side overhead is
#' most of what it costs: swapping one for the other halves the running time.
#' `optim_ignite_adam()` arrived in torch 0.13.0 and DESCRIPTION allows 0.11.0,
#' which is what the fallback is for.
#'
#' @param groups Parameter groups, as `optim_adam()` takes them.
#' @param lr Adam step size.
#' @keywords internal
c2st_adam <- function(groups, lr) {
  if (utils::packageVersion("torch") >= "0.13.0") {
    return(torch::optim_ignite_adam(groups, lr = lr))
  }
  torch::optim_adam(groups, lr = lr)
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
