# c2st() used to shuffle a fold assignment over x and y together
# (base::sample(rep_len(seq_len(n_folds), n))), so at small sample sizes a
# fold's test set could end up with zero rows of one class. roc_auc() returns
# NA_real_ for such a fold, and mean(aucs) (no na.rm) silently made
# c2st()$auc NA with no warning (#269).

test_that("c2st_stratified_folds() gives every fold at least one draw of each class", {
  for (seed in 1:200) {
    set.seed(seed)
    fold <- c2st_stratified_folds(n_x = 6L, n_y = 6L, n_folds = 5L)
    label <- c(rep(0L, 6L), rep(1L, 6L))
    for (k in 1:5) {
      te <- fold == k
      expect_true(any(label[te] == 0L), info = paste("seed", seed, "fold", k))
      expect_true(any(label[te] == 1L), info = paste("seed", seed, "fold", k))
    }
  }
})

test_that("c2st() never returns an NA fold AUC across many seeds at a small sample size", {
  n_each <- 6L
  n_folds <- 5L
  x <- matrix(rnorm(n_each * 2), ncol = 2)
  y <- matrix(rnorm(n_each * 2), ncol = 2)
  for (seed in 1:200) {
    out <- c2st(x, y, n_folds = n_folds, seed = seed, classifier = "logistic")
    expect_false(anyNA(out$fold_auc), info = paste("seed", seed))
    expect_true(is.finite(out$auc), info = paste("seed", seed))
  }
})
